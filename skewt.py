# Copyright (c) 2013-2015 Siphon Contributors.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
"""
Use Siphon to query the NetCDF Subset Service (NCSS).
"""

import datetime as dt

start_time = dt.datetime.now()
import re
from functools import reduce
from metpy.plots import SkewT
import metpy.calc as mpcalc
from metpy.units import units, pandas_dataframe_to_unit_arrays
import matplotlib.pyplot as plt
from siphon.catalog import TDSCatalog
import pandas as pd
import numpy as np
from netCDF4 import num2date

###########################################
# First we construct a TDSCatalog instance using the url
url = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
gfs_catalog = TDSCatalog(url)

# We see this catalog contains three datasets.
# print(gfs_catalog.datasets)
dataset = 'Latest Collection for GFS Quarter Degree Forecast'
gfs_subset = gfs_catalog.datasets[dataset].subset()

###########################################
# Define sub_point to proceed with the query
query = gfs_subset.query()

###########################################
# Then we construct a query asking for data corresponding to desired latitude and longitude and
# for the time interval. We also ask for NetCDF version 4 data and choose the variables.
# This request will return all vertical levels for a single point and for the time interval.
# Note the string representation of the query is a properly encoded query string.
query.lonlat_point(-46, -23)
now = dt.datetime.utcnow()
query.time_range(now, now + dt.timedelta(hours=33))
query.accept('netcdf4')

###########################################
# We'll pull out the variables we want to use, as well as the pressure values.
# To get the name of the correct variable we look at the 'variables' attribute on.
# The last of the variables listed in `coordinates` is the pressure dimension.
# print(gfs_subset.variables)
variable_list = ['Temperature_isobaric',
                 'Relative_humidity_isobaric',
                 'u-component_of_wind_isobaric',
                 'v-component_of_wind_isobaric']
query.variables(*variable_list)

###########################################
# We now request data from the server using this query.
gfs = gfs_subset.get_data(query)


# The variables are then stored in a NetCDF4 dataset
# print(gfs.variables.keys())

def get_variables(data, var_list, time):
    """
    Uses the NetCDF4 dataset from the query,
    the list of variables (with the correct names) and
    time extracted from the variables (gfs.variables['time'][0, :])
    to return a dataset containing the chosen data and pressure levels

    """
    df_list = []
    for variable in var_list:
        chosen_variable = data.variables[variable]
        chosen_variable_vals = chosen_variable[0, time].squeeze()
        press_lvls_name = chosen_variable.coordinates.split()[-1]
        press_lvls = data.variables[press_lvls_name]
        press_lvls_values = press_lvls[0, time].squeeze()
        df = pd.DataFrame([press_lvls_values, chosen_variable_vals], index=None).transpose()
        df.columns = ['pressure', variable]
        df_list.append(df)

    df = reduce(lambda left, right: pd.merge(left, right, on=['pressure'], how='inner'), df_list)
    return df.iloc[::-1]


# Getting time stamps to use in the get_variables() function
times = gfs.variables['time'][0, :]
gfs_output_steps = {}
for time in range(len(times)):
    time_values = gfs.variables['time']
    time_value = num2date(time_values[0, time], time_values.units)
    variables = get_variables(gfs, variable_list, time)
    gfs_output_steps[time_value] = variables # put all df into a dictionary using the timestamp as key


###############
# Adjust the data - change units
def adjust_data(df):
    df['Temperature_isobaric'] = df['Temperature_isobaric'].apply(lambda x: x - 273.15)
    df['pressure'] = df['pressure'].apply(lambda x: x / 100)

    units_dict = {'pressure': 'hPa',
                  'Temperature_isobaric': 'degC',
                  'Relative_humidity_isobaric': 'percent',
                  'u-component_of_wind_isobaric': 'm/s',
                  'v-component_of_wind_isobaric': 'm/s'}

    df = pandas_dataframe_to_unit_arrays(df, units_dict)
    df['Dewpoint'] = mpcalc.dewpoint_from_relative_humidity(df['Temperature_isobaric'],
                                                            df['Relative_humidity_isobaric'])
    return pd.DataFrame(df)


def plot_skewt(df, step):
    """
    takes the treated dataframe and plots the skewt-logp diagram
    """

    # We will pull the data out of the example dataset into individual variables
    # and assign units.
    p = df['pressure'].values * units.hPa
    T = df['Temperature_isobaric'].values * units.degC
    Td = df['Dewpoint'].replace(np.nan, 0.0000001).values * units.degC
    u = df['u-component_of_wind_isobaric'].values * units('meters / second').to('knots')
    v = df['v-component_of_wind_isobaric'].values * units('meters / second').to('knots')

    # Create a new figure. The dimensions here give a good aspect ratio.
    fig = plt.figure(figsize=(12, 9))
    skew = SkewT(fig, rotation=45)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p, T, 'r')
    skew.plot(p, Td, 'g')
    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1013, 100)
    skew.ax.set_xlim(-40, 60)

    # Calculate LCL height and plot as black dot
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

    # Calculate full parcel profile and add to plot as black line
    prof = mpcalc.parcel_profile(p, T[0], Td[0])
    skew.plot(p, prof, 'k', linewidth=2)

    # An example of a slanted line at constant T -- in this case the 0
    # isotherm
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()

    # Add some descriptive titles
    lat = gfs.geospatial_lat_min
    lon = gfs.geospatial_lon_min
    plt.title(f'GFS Sounding ({lat:.2f}, {lon:.2f})', loc='left')

    run_time = re.search(r'(?<=deg_)(.*)(?=.grib2)', gfs.title).group(1)
    plt.title(f'Based on {run_time} GFS run', loc='right')

    valid = step
    plt.title(f'Valid: {valid}', loc='center')

    return skew

# Iterate over the dictionary with all timestamps
# apply the ajust_data() function and generate the skew-t
for step in gfs_output_steps.keys():
    print(step)
    plot_skewt(adjust_data(gfs_output_steps[step]), step)
    plt.savefig(f'./img/{str(step).replace(":", " ").replace(" ", "_")} sounding_gfs.png', format='png')
    plt.show()

end_time = dt.datetime.now() - start_time
print(f'Finished in {end_time.seconds} seconds')

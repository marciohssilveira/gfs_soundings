import re
import numpy as np
from metpy.plots import SkewT
from metpy.units import units
import matplotlib.pyplot as plt


def plot_skewt(adjusted_data, valid):
    """
    :param df: receives the post processed dataframe
    :param valid:
    :return:
    """

    # We will pull the data out of the example dataset into individual variables
    # and assign units.
    p = adjusted_data['pressure'].values * units.hPa
    T = adjusted_data['Temperature_isobaric'].values * units.degC
    Td = adjusted_data['Dewpoint'].replace(np.nan, 0.0000001).values * units.degC
    u = adjusted_data['u-component_of_wind_isobaric'].values * units('meters / second').to('knots')
    v = adjusted_data['v-component_of_wind_isobaric'].values * units('meters / second').to('knots')

    # Create a new figure. The dimensions here give a good aspect ratio.
    fig = plt.figure(figsize=(12, 9))
    skew = SkewT(fig, rotation=45)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p, T, 'r')
    skew.plot(p, Td, 'g')
    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1020, 100)
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

    skew.shade_cape(p, T, prof)
    skew.shade_cin(p, T, prof)

    # Add some descriptive titles
    plt.title(f'GFS Sounding ({station.upper()})', loc='left')

    run_time = re.search(r"(?<=deg_)(.*)(?=.grib2)", data.title).group(1)
    plt.title(f'Based on {run_time} GFS run', loc='right')

    plt.title(f'Valid: {valid}', loc='center')

    return skew

#
# def get_data(station_coordinates, list_of_variables, n_hours):
#     """
#     :param station_coordinates: tuple like (lon, lat)
#     :param list_of_variables: chosen list of variables based on the variables list for the dataset
#     :param n_hours: number of hours for the prediction
#     :return: a subset of the netCDF4 dataset based on the given coordinates and variables
#     """
#     ###########################################
#     # First we construct a TDSCatalog instance using the url
#     URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
#     gfs_catalog = TDSCatalog(URL)
#
#     # We see this catalog contains three datasets.
#     # print(gfs_catalog.datasets)
#     dataset = 'Latest Collection for GFS Quarter Degree Forecast'
#     gfs_subset = gfs_catalog.datasets[dataset].subset()
#
#     ###########################################
#     # Define sub_point to proceed with the query
#     query = gfs_subset.query()
#
#     ###########################################
#     # Then we construct a query asking for data corresponding to desired latitude and longitude and
#     # for the time interval. We also ask for NetCDF version 4 data and choose the variables.
#     # This request will return all vertical levels for a single point and for the time interval.
#     # Note the string representation of the query is a properly encoded query string.
#     query.lonlat_point(*station_coordinates)
#     now = dt.datetime.utcnow()
#     query.time_range(now, now + dt.timedelta(hours=n_hours))
#     query.accept('netcdf4')
#
#     ###########################################
#     # We'll pull out the variables we want to use, as well as the pressure values.
#     # To get the name of the correct variable we look at the 'variables' attribute on.
#     # The last of the variables listed in `coordinates` is the pressure dimension.
#     # print(gfs_subset.variables)
#
#     query.variables(*list_of_variables)
#
#     ###########################################
#     # We now request data from the server using this query.
#     gfs = gfs_subset.get_data(query)
#     return gfs
#
#
# # The variables are then stored in a NetCDF4 dataset
# # print(gfs.variables.keys())
#
# def get_variables(data, list_of_variables, time):
#     """
#     :param data: the subset of netCDF4 obtained with get_data() function
#     :param list_of_variables: chosen list of variables based on the variables list for the dataset
#     :param time:
#     :return:
#     """
#     df_list = []
#     for variable in list_of_variables:
#         chosen_variable = data.variables[variable]
#         chosen_variable_vals = chosen_variable[0, time].squeeze()
#         press_lvls_name = chosen_variable.coordinates.split()[-1]
#         press_lvls = data.variables[press_lvls_name]
#         press_lvls_values = press_lvls[0, time].squeeze()
#         df = pd.DataFrame([press_lvls_values, chosen_variable_vals], index=None).transpose()
#         df.columns = ['pressure', variable]
#         df_list.append(df)
#
#     df = reduce(lambda left, right: pd.merge(left, right, on=['pressure'], how='inner'), df_list)
#     return df.iloc[::-1]
#
#
# ###############
# def adjust_data(df):
#     """
#     :param df: receives the dataframe containing the data and adjust the units
#     :return: the updated dataframe
#     """
#     # converts temperature from Kelvin to Degrees Celsius
#     df['Temperature_isobaric'] = df['Temperature_isobaric'].apply(lambda x: x - 273.15)
#
#     # converts pressure from Pa to hPa
#     df['pressure'] = df['pressure'].apply(lambda x: x / 100)
#
#     # dictionary with the units used to apply into the mpcalc function
#     # the package metpy functions will only deal with data with its units associated
#     units_dict = {'pressure': 'hPa',
#                   'Temperature_isobaric': 'degC',
#                   'Relative_humidity_isobaric': 'percent',
#                   'u-component_of_wind_isobaric': 'm/s',
#                   'v-component_of_wind_isobaric': 'm/s'}
#
#     # Attach units to data into the dataframe and return united arrays
#     df = pandas_dataframe_to_unit_arrays(df, units_dict)
#
#     # Calculate the ambient dewpoint given air temperature and relative humidity.
#     df['Dewpoint'] = mpcalc.dewpoint_from_relative_humidity(df['Temperature_isobaric'],
#                                                             df['Relative_humidity_isobaric'])
#
#     # converto to pandas dataframe again as the plt_skew() metpy function suggests
#     return pd.DataFrame(df)
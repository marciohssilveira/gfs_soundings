"""
Skew-T Analysis
===============

Classic skew-T/log-p plot using data from University of Wyoming.

This example uses example data from the University of Wyoming sounding
archive for 12 UTC 31 October 2016 for Minneapolis, MN (MPX) and uses
MetPy to plot the classic skew-T with Temperature, Dewpoint, and wind
barbs.

"""

from datetime import datetime

import matplotlib.pyplot as plt
from metpy.plots import SkewT
from metpy.units import pandas_dataframe_to_unit_arrays, units
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir


######################################################################
# Set time using a datetime object and station as variables
#

dt = datetime(2020, 7, 13, 12)
station = 'SBMT'


######################################################################
# Grab Remote Data
# ----------------
#
# This requires an internet connection to access the sounding data from a
# remote server at the University of Wyoming.
#

# Read remote sounding data based on time (dt) and station
data = WyomingUpperAir.request_data(dt, station)


def plot_skewt(df):
    # We will pull the data out of the example dataset into individual variables
    # and assign units.
    p = df['pressure'].values * units.hPa
    T = df['temperature'].values * units.degC
    Td = df['dewpoint'].replace(np.nan, 0.0000001).values * units.degC
    u = df['u_wind'].values * units('meters / second')
    v = df['v_wind'].values * units('meters / second')

    # Create a new figure. The dimensions here give a good aspect ratio.
    fig = plt.figure(figsize=(11, 9))
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
    plt.title(f'{station} Sounding', loc='left')
    plt.title(f'Valid Time: {dt}', loc='right')

    return skew


plot_skewt(data)

plt.savefig('sounding_real.png', format='png')
plt.show()


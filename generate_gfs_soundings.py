import datetime as dt

start = dt.datetime.now()
from get_gfs_data import GetGFSData
from plot_skew import PlotSkew
from netCDF4 import num2date
import matplotlib.pyplot as plt
import re
import os
import warnings

warnings.filterwarnings("ignore")

coordinates = {'SBGR': (-46.473056, -23.435556),
               'SBMT': (-46.633889, -23.506667),
               'SBSP': (-46.656389, -23.626111),
               'SBJD': (-46.943611, -23.181667),
               'SBBP': (-46.537500, -22.979167),
               'SBKP': (-47.134444, -23.006944),
               'SDCO': (-47.486389, -23.483056),
               'SBJH': (-47.165833, -23.426944),
               'SBAQ': (-48.140278, -21.804444),
               'SBRP': (-47.776667, -21.136389),
               'SBSR': (-49.404722, -20.816111),
               'SBUR': (-47.966111, -19.764722),
               'SBUL': (-48.225278, -18.883611),
               'SBAX': (-46.965556, -19.560556),
               'SBVG': (-45.473333, -21.588889)}

n_hours = 3
URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset = 'Latest Collection for GFS Quarter Degree Forecast'
variables = ['Temperature_isobaric',
             'Relative_humidity_isobaric',
             'u-component_of_wind_isobaric',
             'v-component_of_wind_isobaric']
# pressure levels are extracted from each variable
units = ['hPa', 'degC', 'percent', 'm/s', 'm/s']

skewt = GetGFSData(URL, dataset, variables, units)

folders_to_be_made = list(coordinates.keys())
root_path = './img'
for folder in folders_to_be_made:
    if folder not in os.listdir(root_path):
        os.mkdir(os.path.join(root_path, folder))

for station, coordinate in coordinates.items():
    raw_data = skewt.get_gfs_data(station, coordinate, n_hours)
    time_steps = raw_data.variables['time'][0, :]
    for time_step in range(len(time_steps)):
        step_values = raw_data.variables['time']
        step_value = num2date(step_values[0, time_step], step_values.units)
        station_data = skewt.get_variables(raw_data, time_step)
        PlotSkew(station_data).plot_skewt()
        # Add some descriptive titles
        plt.title(f'GFS Sounding ({station.upper()})', loc='left')
        run_time = re.findall(r'(?<=deg_)(.*)(?=.grib2)', raw_data.title)[0]
        plt.title(f'Based on {run_time} GFS run', loc='right')
        plt.title(f'Valid: {step_value}', loc='center')
        step = f'{step_value.day:02}_{step_value.hour:02}'
        plt.savefig(f'./img/{station.upper()}/{station.upper()}_GFS_sounding_{step}UTC.png', format='png')

end = dt.datetime.now()
print('Done!')
print(f'Total time elapsed: {end - start}')

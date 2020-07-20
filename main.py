from get_data import GetGFSData
from netCDF4 import num2date

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

n_hours = 24
URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset = 'Latest Collection for GFS Quarter Degree Forecast'
variables = ['Temperature_isobaric',
                 'Relative_humidity_isobaric',
                 'u-component_of_wind_isobaric',
                 'v-component_of_wind_isobaric']
# pressure levels are extracted from each variable
units = ['hPa', 'degC', 'percent', 'm/s', 'm/s']


skewt = GetGFSData(URL, dataset, variables, units)


for station, coordinate in coordinates.items():
    raw_data = skewt.get_gfs_data(coordinate, n_hours)
    time_steps = raw_data.variables['time'][0, :]
    raw_data_time_steps = {}
    for time_step in range(len(time_steps)):
        step_values = raw_data.variables['time']
        step_value = num2date(step_values[0, time_step], step_values.units)
        variables = skewt.get_variables(raw_data, time_step)
        raw_data_time_steps[step_value] = variables  # put all df into a dictionary using the timestamp as key

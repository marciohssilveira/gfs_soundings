import warnings
import os
import re
import matplotlib.pyplot as plt
from netCDF4 import num2date
from skewt_logp import PlotSkew
from get_gfs_data import GetGFSData
import datetime as dt

start = dt.datetime.now()

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    def plot_skewt(coordinates, n_hours, gfs_data):
        for station, coordinate in coordinates.items():
            raw_data = gfs_data.get(station, coordinate, n_hours)
            time_steps = raw_data.variables['time'][0, :]
            for time_step in range(len(time_steps)):
                step_values = raw_data.variables['time']
                step_value = num2date(
                    step_values[0, time_step], step_values.units)
                station_data = gfs_data.get_variables(raw_data, time_step)
                PlotSkew().plot_skewt(station_data)
                # Add some descriptive titles
                plt.title(f'GFS Sounding ({station.upper()})', loc='left')
                run_time = re.findall(
                    r'(?<=deg_)(.*)(?=.grib2)', raw_data.title)[0]
                plt.title(f'Based on {run_time} GFS run', loc='right')
                plt.title(f'Valid: {step_value}', loc='center')
                step = f'{step_value.day:02}_{step_value.hour:02}'
                plt.savefig(
                    f'./img/{station.upper()}_GFS_sounding_{step}UTC.png', format='png')

        end = dt.datetime.now()
        print('Done!')
        print(f'Total time elapsed: {end - start}')


# defining coordinates and other variables

    area_1 = {'SBSR': (-49.404722, -20.816111),
              'SBAU': (-50.426389, -21.144167),
              'SBAQ': (-48.140278, -21.804444),
              'SBRP': (-47.776667, -21.136389),
              'SBAX': (-46.965556, -19.560556),
              'SBUL': (-48.225278, -18.883611),
              'SBUR': (-47.966111, -19.764722),
              'SBVG': (-45.473333, -21.588889),
              'SBZM': (-43.173056, -21.513056),
              'SBBH': (-43.950556, -19.851944),
              'SBCF': (-43.971944, -19.624444),
              'SBMK': (-43.821944, -16.706111),
              'SBIP': (-42.488056, -19.470556),
              'SBGV': (-41.986111, -18.896944),
              'SBBR': (-47.918611, -15.871111),
              'SBGO': (-49.221111, -16.632500),
              'SBCN': (-48.610000, -17.724722),
              'SWLC': (-50.956111, -17.834722),
              'SBBW': (-52.389444, -15.860833),
              'SNBR': (-45.009444, -12.079167),
              'SBLE': (-41.276944, -12.482222),
              'SBVC': (-40.914722, -14.907778),
              'SNVB': (-38.992500, -13.296389),
              'SDIY': (-38.900556, -12.200556),
              'SBSV': (-38.322500, -12.908611),
              'SBIL': (-39.033333, -14.815000),
              'SNTF': (-39.668333, -17.524444)}

    area_2 = {'SBGR': (-46.473056, -23.435556),
              'SBMT': (-46.633889, -23.506667),
              'SBSP': (-46.656389, -23.626111),
              'SBSJ': (-45.871111, -23.228889),
              'SBTA': (-45.515833, -23.038889),
              'SBST': (-46.299722, -23.928056),
              'SBKP': (-47.134444, -23.006944),
              'SBJD': (-46.943611, -23.181667),
              'SBBP': (-46.537500, -22.979167),
              'SBJH': (-47.165833, -23.426944),
              'SDCO': (-47.486389, -23.483056),
              'SBBU': (-49.053889, -22.343611),
              'SBAE': (-49.068333, -22.157778),
              'SBML': (-49.926944, -22.195556),
              'SBDN': (-51.418889, -22.178333),
              'SBCR': (-57.671389, -19.011944),
              'SBCG': (-54.670278, -20.469444),
              'SBTG': (-51.680278, -20.751389),
              'SBDB': (-56.452500, -21.247222),
              'SBDO': (-54.925556, -22.200556),
              'SBPP': (-55.703056, -22.549722),
              'SBLO': (-51.136667, -23.330278),
              'SBMG': (-52.012222, -23.479444),
              'SBTD': (-53.696389, -24.685278),
              'SBCA': (-53.501944, -25.002222),
              'SSGG': (-51.523611, -25.388333),
              'SBPO': (-52.694444, -26.217222)}

    n_hours = 3

    gfs_data = GetGFSData()

    plot_skewt(area_1, n_hours, gfs_data)
    plot_skewt(area_2, n_hours, gfs_data)

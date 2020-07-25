from cx_Freeze import setup, Executable

setup(name = 'gfs_soundings',
      version = '1.0' ,
      description = "Plotting Skew-T/Log-P diagrams using GFS output data" ,
      executables = [Executable("generate_gfs_soundings.py")])
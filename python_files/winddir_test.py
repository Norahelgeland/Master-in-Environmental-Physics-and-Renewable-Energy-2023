import numpy as np
import pandas as pd
from datetime import datetime
import time
import csv
import readline
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import os
#from netCDF4 import Dataset, num2date, date2num
import math
import sys
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import xarray
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import gc
from numpy import datetime64
import fimex_xarray
from netCDF4 import Dataset
import time
from radiosonde_class import Radiosonde
from tqdm import tqdm
import os
import matplotlib
matplotlib.use('QtAgg')



def get_projection(ds):

    if "proj4" in ds.attrs:
        proj_string = ds.attrs["proj4"].removeprefix("+init=")
        if proj_string == "epsg:32633":
            return ccrs.UTM(33)
        # return ccrs.UTM(33)
        if proj_string == "+proj=longlat":
            return ccrs.PlateCarree()
        return ccrs.Projection(proj_string)
    pr = ds["projection_lambert"].attrs
    if pr["grid_mapping_name"] == "rotated_latitude_longitude":
        projection = ccrs.RotatedPole(
            pole_longitude=pr["grid_north_pole_longitude"],
            pole_latitude=pr["grid_north_pole_latitude"],
        )
    elif pr["grid_mapping_name"] == "latitude_longitude":
        projection = ccrs.PlateCarree()
    elif pr["grid_mapping_name"] == "lambert_conformal_conic":
        standard_parallel = pr["standard_parallel"]
        try:
            standard_parallel[1]
        except IndexError:
            standard_parallel = [standard_parallel, standard_parallel]
        projection = ccrs.LambertConformal(
            central_longitude=pr["longitude_of_central_meridian"],
            central_latitude=pr["latitude_of_projection_origin"],
            standard_parallels=standard_parallel,
        )
    else:
        raise NotImplementedError(f"Projection {pr['grid_mapping_name']} not known")

    return projection


if __name__=="__main__":

    config = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/netcdf/cdmGribReaderConfigArome2_5.xml"
    arome_file = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/gribml/deter_19860425_00.grbml"

    ds = xarray.open_dataset(arome_file, config=config, engine="fimex")

    ds = ds.isel(x=slice(None, None, 100), y=slice(None, None, 100))
    projection=get_projection(ds)

    fig, ax = plt.subplots(1,subplot_kw={"projection": projection})
    #ax.pcolormesh(ds["x"], ds["y"], ds["air_temperature_2m"].isel(time=1).squeeze(), transform=projection)
    ax.quiver(ds["x"], ds["y"], ds["x_wind_ml"].isel(time=1, hybrid=-1), ds["y_wind_ml"].isel(time=1, hybrid=-1), scale = 10, transform=projection)
    ax.coastlines()
    ax.set_extent((-10, 60, 40, 75), crs=ccrs.PlateCarree())
   # ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.LAKES)
    #ax.add_feature(cfeature.RIVERS)
    #ax.add_feature(cfeature.BORDERS)

    #ax.gridlines(draw_labels=["left", "bottom"], x_inline=False, y_inline=False)
    file = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/alpha_chernobyl.nc"
    ds_alpha = xarray.open_dataset(file)   
    ds_alpha = ds_alpha["alpha"]
    print("x_wind")
    print(ds["x_wind_ml"].isel(hybrid=63, x=5, y=5, time=1).values)
    print("y_wind_ml")
    print(ds["y_wind_ml"].isel(hybrid=63, x=5, y=5, time=1).values)
    print("alpha")
    print(ds_alpha.isel(x=5, y=5).values)

    ds = ds.isel(time=1)
    x, y = np.meshgrid(ds["x"], ds["y"])

    #stuff = projection.transform_vectors(ccrs.PlateCarree(), np.ravel(x), np.ravel(y), np.ravel(ds["x_wind_ml"].isel(hybrid=-1)), np.ravel(ds["y_wind_ml"].isel(hybrid=-1)))
    lonlat_projection = ccrs.PlateCarree()
    
    stuff = lonlat_projection.transform_points(projection, np.ravel(x), np.ravel(y))    
    lon = stuff[:, 0]
    lat = stuff[:, 1]
    vlon, vlat = lonlat_projection.transform_vectors(projection, np.ravel(x), np.ravel(y), np.ravel(ds["x_wind_ml"].isel(hybrid=-1)), np.ravel(ds["y_wind_ml"].isel(hybrid=-1)))

    ax.quiver(lon, lat, vlon, vlat, color="b", transform=lonlat_projection, scale=10)

    #breakpoint()
   # wdir = ds_alpha+90-np.arctan2(ds["y_wind_ml"].isel(time=i).values, ds["x_wind_ml"].isel(time=1).values)*180.0/math.pi
 
    plt.show()
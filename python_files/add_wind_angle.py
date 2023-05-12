
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

    # loading the Chernobyl data from the Arome model
    config = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/netcdf/cdmGribReaderConfigArome2_5.xml"
    arome_file = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/gribml/deter_19860425_00.grbml"

    ds = xarray.open_dataset(arome_file, config=config, engine="fimex")
    file = "/lustre/storeB/project/fou/kl/cerad/Meteorology/AROME-CHERNOBYL/alpha_chernobyl.nc"

    ds_alpha = xarray.open_dataset(file)   
    ds_alpha = ds_alpha["alpha"]
    projection=get_projection(ds)   
    #alpha_list = [ds_alpha]*19
    #ds_alpha = xarray.concat(alpha_list, "time")

    #alpha_list2 = [ds_alpha]*65
    #ds_alpha = xarray.concat(alpha_list2, "hybrid")

    i = int(os.getenv("SGE_TASK_ID")) - 1
    if True:

        ds2 = ds.isel(time=i)
        x, y = np.meshgrid(ds2["x"], ds2["y"])

        #stuff = projection.transform_vectors(ccrs.PlateCarree(), np.ravel(x), np.ravel(y), np.ravel(ds["x_wind_ml"].isel(hybrid=-1)), np.ravel(ds["y_wind_ml"].isel(hybrid=-1)))
        lonlat_projection = ccrs.PlateCarree()
    
        stuff = lonlat_projection.transform_points(projection, np.ravel(x), np.ravel(y))    
        lon = stuff[:, 0]
        lat = stuff[:, 1]

        vlon, vlat = lonlat_projection.transform_vectors(projection, np.ravel(x), np.ravel(y), np.ravel(ds2["x_wind_ml"].isel(hybrid=0)), np.ravel(ds2["y_wind_ml"].isel(hybrid=0)))
        vlon = np.reshape(vlon, x.shape)
        vlat = np.reshape(vlat, y.shape)
        #wdir = ds_alpha+90-np.arctan2(ds["y_wind_ml"].isel(time=i).values, ds["x_wind_ml"].isel(time=i).values)*180.0/math.pi
        wdir = 90-np.arctan2(vlon, vlat)*180.0/math.pi

        # angles=np.arctan2(ds["y_wind_ml"].isel(time=i).values, ds["x_wind_ml"].isel(time=i).values)*180.0/math.pi

        wind_dir_data = np.zeros((65, len(ds["y"]), len(ds["x"])))

        wind_dir_data[0, :, :] = wdir

        for h in range(1,65):

            ds3 = ds2.copy()
            x, y = np.meshgrid(ds3["x"], ds3["y"])

            #stuff = projection.transform_vectors(ccrs.PlateCarree(), np.ravel(x), np.ravel(y), np.ravel(ds["x_wind_ml"].isel(hybrid=-1)), np.ravel(ds["y_wind_ml"].isel(hybrid=-1)))
            lonlat_projection = ccrs.PlateCarree()
    
            stuff = lonlat_projection.transform_points(projection, np.ravel(x), np.ravel(y))    
            lon = stuff[:, 0]
            lat = stuff[:, 1]
            vlon, vlat = lonlat_projection.transform_vectors(projection, np.ravel(x), np.ravel(y), np.ravel(ds3["x_wind_ml"].isel(hybrid=h)), np.ravel(ds3["y_wind_ml"].isel(hybrid=h)))

            #wdir = ds_alpha+90-np.arctan2(ds["y_wind_ml"].isel(time=i).values, ds["x_wind_ml"].isel(time=i).values)*180.0/math.pi
            vlon = np.reshape(vlon, x.shape)
            vlat = np.reshape(vlat, y.shape)
            wdir = 90-np.arctan2(vlon, vlat)*180.0/math.pi

            wind_dir_data[h, :, :] = wdir

            #windir_data2=xarray.DataArray(wdir, dims=("y", "x"))

            # angles=np.arctan2(ds["y_wind_ml"].isel(time=i).values, ds["x_wind_ml"].isel(time=i).values)*180.0/math.pi
        
            #windir_data2["wdir"] = wdir
     
        windir_data = xarray.DataArray(wind_dir_data, dims=["hybrid", "y", "x"], coords={"x": ds3["x"], "y":ds3["y"], "hybrid":ds["hybrid"]})

        windir_data = xarray.Dataset({"wdir": windir_data})
        windir_data.to_netcdf("windir_data_"+str(i)+".nc")


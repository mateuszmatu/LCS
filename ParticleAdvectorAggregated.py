import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from datetime import timedelta, datetime
import xarray as xr
import matplotlib.pyplot as plt
import os



class AggregatedBackwardsAdvection:
    """
        Class for advecting a grid of particles backwards in time using the aggregated Barents-2.5 EPS file from 
        'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
    """
    def __init__(self, lons, lats, ts, sep, dur, date, at_time=0):
        """
            Initiates the class. 
        Args:
            lons    [list]     :      List of two lon values, [lon1, lon2], where lon1 is in the bottom left corner of domain and lon2 is in the top right corner of domain.
            lats    [list]     :      List of two lat values, [lat1, lat2], where lat1 is in the bottom left corner of domain and lat2 is in the top right corner of domain.
            ts      [float]    :      The integration time step. Set negative for forwards in time advection.
            sep     [int]      :      Initial separation between gridded particles in meters. 
            date    [list]     :      Date on which particles are seeded. [year, day, month].
            at_time [int]      :      Hour of the day for the particles to be initiated on.  
        """
        self.lons=lons
        self.lats=lats
        #self.proj = '+proj=lcc +lat_0=77.5 +lon_0=-25 +lat_1=77.5 +lat_2=77.5 +no_defs +R=6.371e+06'
        self.ts=ts
        self.sep=sep
        self.dur=dur
        self.date = datetime(date[0], date[1], date[2], at_time)
        #self.file = 'Means/daily_mean_velocity_april_m1.nc'
        self.file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        self.at_time=at_time

    def advect(self, outfile):
        """
            Advects grid of particles backwards in time (forwards in time if ts < 0) using the OpenDrift module developed by MET Norway. 
            Saves the initial grid and final positions of particles to file.
        Args:
            outfile     [str]   :   name of output file. 
        """
        o = OceanDrift(loglevel=20)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        r = reader_netCDF_CF_generic.Reader(self.file) #proj4=self.proj)

        x, y = r.lonlat2xy(self.lons, self.lats)
        c1 = np.arange(x[0], x[1], self.sep)
        c2 = np.arange(y[0], y[1], self.sep)
        X, Y = np.meshgrid(c1, c2)

        lons, lats = r.xy2lonlat(X.flatten(), Y.flatten())

        o.add_reader(r)
        o.seed_elements(lons.ravel(), lats.ravel(), time=self.date, radius_type='uniform')
        o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=-self.ts), time_step_output=timedelta(hours=self.dur), outfile=f'{outfile}.nc')

        lons, lats = np.reshape(lons, (X.shape[0], X.shape[1])), np.reshape(lats, (X.shape[0], X.shape[1]))
        f_x1, f_y1 = r.lonlat2xy(o.history['lon'].T[-1], o.history['lat'].T[-1])

        d = xr.open_dataset(f'{outfile}.nc')

        os.remove(f'{outfile}.nc')
        ds = xr.Dataset(coords=dict(lon = (['x', 'y'], X),
                        lat = (['x','y'], Y)),
                data_vars=dict(separation=self.sep, duration=self.dur, nlon=f_x1, nlat=f_y1))
        ds.to_netcdf(f'{outfile}.nc')



class AdvectForwardsGrid:
    """
        Class for advecting a grid of particles forwards in time using the aggregated Barents-2.5 EPS file from 
        'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
    """
    def __init__(self, lons, lats, ts, sep, dur, date, at_time=0):
        """
            Initiates the class. 
        Args:
            lons    [list]     :      List of two lon values, [lon1, lon2], where lon1 is in the bottom left corner of domain and lon2 is in the top right corner of domain.
            lats    [list]     :      List of two lat values, [lat1, lat2], where lat1 is in the bottom left corner of domain and lat2 is in the top right corner of domain.
            ts      [float]    :      The integration time step. Set negative for backwards in time advection.
            sep     [int]      :      Initial separation between gridded particles in meters. 
            date    [list]     :      Date on which particles are seeded. [year, day, month].
            at_time [int]      :      Hour of the day for the particles to be initiated on.  
        """
        self.lons=lons
        self.lats=lats
        self.proj = '+proj=lcc +lat_0=77.5 +lon_0=-25 +lat_1=77.5 +lat_2=77.5 +no_defs +R=6.371e+06'
        self.ts=ts
        self.sep=sep
        self.dur=dur
        self.date = datetime(date[0], date[1], date[2], at_time)
        self.file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
    
    def advect(self, outfile):
        """
            Advects a grid of particles forwards in time, and saves all intermediate positions of particles to file. 
        Args:
            outfile     [str]   :   name of output file
        """
        o = OceanDrift(loglevel=20)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        r = reader_netCDF_CF_generic.Reader(self.file, proj4=self.proj)

        x, y = r.lonlat2xy(self.lons, self.lats)
        c1 = np.arange(x[0], x[1], self.sep)
        c2 = np.arange(y[0], y[1], self.sep)
        X, Y = np.meshgrid(c1, c2)

        lons, lats = r.xy2lonlat(X.flatten(), Y.flatten())

        o.add_reader(r)

        o.seed_elements(lons.ravel(), lats.ravel(), time=self.date, radius_type='uniform')
        o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=self.ts), outfile=f'{outfile}.nc')


test = AdvectForwardsGrid(lons=[10, 14], lats=[67.3,67.7], ts = 3600, dur=1, number=25000, sep = 500, radius=30000, date=[2023,2,1], at_time=12)
test.advect('particle_test_grid_2')






        
        




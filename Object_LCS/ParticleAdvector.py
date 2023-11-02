
import functions as f
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from datetime import timedelta, datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import os
from opendrift.readers import reader_double_gyre


class AdvectionBarentsEPS:
    """
        Class for advecting a grid of particles in an ensemble member from the Barents-2.5 EPS. Internal use only. 
    """
    def __init__(self, lons, lats, ts, sep, dur, date, at_time=0):
        """
            Initiates the class. 
        Args:
            lons    [list]     :      List of two lon values, [lon1, lon2], where lon1 is in the bottom left corner of domain and lon2 is in the top right corner of domain.
            lats    [list]     :      List of two lat values, [lat1, lat2], where lat1 is in the bottom left corner of domain and lat2 is in the top right corner of domain.
            ts      [float]    :      The integration time step. Set negative for backwards in time advection. In Seconds
            sep     [int]      :      Initial separation between gridded particles in meters. 
            date    [str]      :      Date on which particles are seeded. 'ymd' without spacings.
            at_time [int]      :      Hour of the day for the particles to be initiated on.  
        """
        self.lons=lons
        self.lats=lats
        self.proj = '+proj=lcc +lat_0=77.5 +lon_0=-25 +lat_1=77.5 +lat_2=77.5 +no_defs +R=6.371e+06'
        self.ts=ts
        self.sep=sep
        self.dur=dur
        self.date=date
        self.at_time=at_time

        self.tf = f.files_from_lustre(date)
        self.name = f.name_from_lon_lat(lons, lats)

    def displace_one_member(self, member, outfile=None):
        """
            Advects gridded particles using velocity field from one ensemble member. Save initiall and final positions of particles to file
        Args:
            member     [int]    :   The ensemble member from which velocity field is to be used. Number between 0-23.
            outfile    [str]    :   Name of file output is saved to.
        """
        ha = f.hour_adjustment(member) + self.at_time
        corr = f.correct_file(self.tf, member)
        file = corr[0]
        _m = corr[1]

        #Defining the savefile location
        if outfile is None:
            string = f'{self.path}/{self.name}_h{self.at_time}_m{member}'
        else:
            string = f'{outfile}'

        o = OceanDrift(loglevel=20)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_m, proj4=self.proj)

        x, y = r.lonlat2xy(self.lons, self.lats)
        c1 = np.arange(x[0], x[1], self.sep)
        c2 = np.arange(y[0], y[1], self.sep)
        X, Y = np.meshgrid(c1, c2)

        lons, lats = r.xy2lonlat(X.flatten(), Y.flatten())
        
        #Starting OpenDrift
        
        o.add_reader(r)

        o.seed_elements(lons.ravel(), lats.ravel(), time=r.start_time+timedelta(hours=ha), radius_type='uniform')
        o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=self.ts), time_step_output=timedelta(hours=self.dur), outfile=f'{string}.nc')
        lons, lats = np.reshape(lons, (X.shape[0], X.shape[1])), np.reshape(lats, (X.shape[0], X.shape[1]))
        f_x1, f_y1 = r.lonlat2xy(o.history['lon'].T[-1], o.history['lat'].T[-1])
        
        d = xr.open_dataset(f'{string}.nc')

        os.remove(f'{string}.nc')
        ds = xr.Dataset(coords=dict(lon = (['x', 'y'], X),
                        lat = (['x','y'], Y)),
                data_vars=dict(separation=self.sep, duration=self.dur, nlon=f_x1, nlat=f_y1))
        ds.to_netcdf(f'{string}.nc')

        return string

class Advection:
    def __init__(self, file, lons, lats, ts, sep, dur, start):
        '''
        Args:
            DataArray   (xr.DataArray)  :   A file containing the model data.
            lons    (list)              :   A list containing max and min lon
            lats    (list)              :   A list containing max and min lat
            ts      (int)               :   Simulation timestep
            sep     (int)               :   initial separation between grid particles in meters.
            dur     (int)               :   simulation duration in hours.
            start   (Datetime)          :   A datetime object containing the simulation start time.
        '''
        self.file = file
        self.lons = lons
        self.lats = lats
        self.ts = ts
        self.sep = sep
        self.dur = dur
        self.start = start
    
    def run(self, outfile = None, ensemble_member=None):
        o = OceanDrift(loglevel=20)
        try:
            r = reader_netCDF_CF_generic.Reader(self.file, ensemble_member=ensemble_member)
        except:
            pass
            #r = reader_ROMS_native.Reader(self.file)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        o.set_config('drift:vertical_mixing', False)
        o.add_reader(r)

        x,y = r.lonlat2xy(self.lons, self.lats)
        if x[1]>x[0]:
            c1 = np.arange(x[0], x[1], self.sep)
        else:
            c1 = np.arange(x[0], x[1], -self.sep)
        if y[1]>y[0]:
            c2 = np.arange(y[0], y[1], self.sep)
        else:
            c2 = np.arange(y[0], y[1], -self.sep)
        X, Y = np.meshgrid(c1, c2)
  
        lons, lats = r.xy2lonlat(X.flatten(), Y.flatten())

        o.seed_elements(lons.ravel(), lats.ravel(), time=r.start_time+self.start)
        o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=self.ts), time_step_output=timedelta(hours=self.dur))
        lons, lats = np.reshape(lons, (X.shape[0], X.shape[1])), np.reshape(lats, (X.shape[0], X.shape[1]))
        f_x1, f_y1 = r.lonlat2xy(o.history['lon'].T[-1], o.history['lat'].T[-1])


        ds = xr.Dataset(coords=dict(lon = (['x', 'y'], X),
                        lat = (['x','y'], Y)),
                data_vars=dict(separation=self.sep, duration=self.dur, nlon=f_x1, nlat=f_y1))
        
        if outfile is not None:
            ds.to_netcdf(f'{outfile}.nc')

        return ds
        
class DoubleGyre:
    """
        Class for advecting a grid of particles using velocity fields from the analytical double gyre system .
    """
    def __init__(self, at_time = 3, dur = 15, time_step = 0.5, sep = 0.01, epsilon = 0.25, omega = 0.682, A = 0.1, outfile = None):
        """
            Initiates class
        Args:
            at_time     [float]     :   Time for which particles should be initiated. 
            dur         [int]       :   Duration of particle advection.
            time_step   [flaot]     :   Time step for the integration. Set negative for backwards in time integration.
            sep         [float]     :   Separation between particles. Smaller numbers means smaller separation between particles. 
            epsilon     [float]     :   The epsilon parameter of the double gyre. Expansion and contraction. 
            omega       [float]     :   Omega parameter. Frequency of oscillations. 
            A           [float]     :   Rotational velocity of the two gyres.
            outfile     [str]       :   File where output is saved to.
        """
        self.at = at_time
        self.dur = dur
        self.ts = time_step
        self.sep = sep
        self.eps = epsilon
        self.om = omega
        self.A = A 
        self.outfile=outfile
    
    def advect(self):
        """
            Advects the grid of particles forwards/backwards in time (based on sign of ts). 
            Saves initial and final positions of the gridded particles to file. 
        """
        o = OceanDrift(loglevel=30)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        double_gyre = reader_double_gyre.Reader(epsilon=self.eps, omega=self.om, A=self.A)
        o.add_reader(double_gyre)

        proj = double_gyre.proj
        x = np.arange(double_gyre.xmin, double_gyre.xmax, self.sep)
        y = np.arange(double_gyre.ymin, double_gyre.ymax, self.sep)

        X, Y = np.meshgrid(x,y)
        lons, lats = proj(X,Y, inverse=True)

        o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=self.at))
        o.run(duration=timedelta(seconds=self.dur), time_step=self.ts, time_step_output=self.dur, outfile=f'{self.outfile}.nc')
        f_x1, f_y1 = proj(o.history['lon'].T[-1], o.history['lat'].T[-1])
        string = self.outfile

        d = xr.open_dataset(f'{string}.nc')
        os.remove(f'{string}.nc')
        ds = xr.Dataset(coords=dict(lon = (['x', 'y'], X),
                                lat = (['x','y'], Y)),
                        data_vars=dict(separation=self.sep, duration=self.dur, nlon=f_x1, nlat=f_y1))

        ds.to_netcdf(f'{string}.nc')




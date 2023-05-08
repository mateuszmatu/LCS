import functions as f
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from datetime import timedelta, datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import os
from opendrift.readers import reader_double_gyre


file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps/barents_eps_20221125T00Z.nc'

def perturbed_position():
    """
        Advects two initially close particles in a velocity field. 
        Saves their trajectories to file in 'random_files/perturbed_position.nc' and plots the trajectories. 
    """
    o = OceanDrift(loglevel=20)
    o.set_config('drift:advection_scheme', 'runge-kutta4')
    r = reader_netCDF_CF_generic.Reader(file, ensemble_member=0)
    o.add_reader(r)
    lon = [15.6, 15.8]
    lat = [70.2, 70.2]

    o.seed_elements(lon, lat, time=r.start_time, number=2)
    o.run(duration=timedelta(hours=48), time_step=3600, outfile='random_files/perturbed_position.nc')
    o.plot()

def perturbed_velocity():
    """
        Advects two particles initially positioned at the same location with different velocity fields. 
        Saves their trajectories to file in 'random_files/perturbed_velocity{i}.nc', where i=0,1, and plots the trajectories. 
    """
    for i in range(2):
        o = OceanDrift(loglevel=20)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        r = reader_netCDF_CF_generic.Reader(file, ensemble_member=i+3)
        o.add_reader(r)
        lon = 15.4
        lat = 70.2
        
        o.seed_elements(lon, lat, time=r.start_time)
        o.run(duration=timedelta(hours=48), time_step=3600, outfile=f'random_files/perturbed_velocity{i}.nc')
    
    pos1 = xr.open_dataset('random_files/perturbed_velocity0.nc')
    pos2 = xr.open_dataset('random_files/perturbed_velocity1.nc')

    plt.scatter(pos1.lon, pos1.lat, color='blue')
    plt.scatter(pos2.lon, pos2.lat, color='red')
    plt.show()

def both_perturbed():
    """
        Advects two particles initially positioned some distance away from each other using two different velocity fields. 
        Saves their trajectories to file in 'random_files/perturbed_both{i}.nc', where i=0,1, and plots the trajectories. 
    """
    lon = [15.6, 15.8]
    lat = [70.2, 70.2]
    for i in range(2):
        o = OceanDrift(loglevel=20)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        r = reader_netCDF_CF_generic.Reader(file, ensemble_member=i+2)
        o.add_reader(r)
        
        o.seed_elements(lon[i], lat[i], time=r.start_time)
        o.run(duration=timedelta(hours=48), time_step=3600, outfile=f'random_files/perturbed_both{i}.nc')
    
    pos1 = xr.open_dataset('random_files/perturbed_both0.nc')
    pos2 = xr.open_dataset('random_files/perturbed_both1.nc')

    plt.scatter(pos1.lon, pos1.lat, color='blue')
    plt.scatter(pos2.lon, pos2.lat, color='red')
    plt.show()

def plot_all():
    """
        Reads in files from the three functions above, and plots the trajectories of all particles in three separate subplots. 
        Saves the figure to 'figures/particle_uncertainty.pdf.
    """
    pos1 = xr.open_dataset('random_files/perturbed_position.nc')
    pos2 = xr.open_dataset('random_files/perturbed_velocity0.nc')
    pos3 = xr.open_dataset('random_files/perturbed_velocity1.nc')
    pos4 = xr.open_dataset('random_files/perturbed_both0.nc')
    pos5 = xr.open_dataset('random_files/perturbed_both1.nc')
    fig, axs = plt.subplots(1,3,figsize=(10,4.5),subplot_kw={'projection':ccrs.NorthPolarStereo(central_longitude=-24.05)})
    #fig.set_aspect('equal')
    axs[0].scatter(pos1.lon[0,:], pos1.lat[0,:], color='black', transform=ccrs.PlateCarree())
    axs[0].scatter(pos1.lon[1,:], pos1.lat[1,:], color='black', transform=ccrs.PlateCarree())

    axs[1].scatter(pos2.lon, pos2.lat, color='black', transform=ccrs.PlateCarree())
    axs[1].scatter(pos3.lon, pos3.lat, color='black', transform=ccrs.PlateCarree())
    #axs[1].set_aspect('equal', crs=ccrs.NorthPolarStereo(central_longitude=-24.05))
    axs[2].scatter(pos4.lon, pos4.lat, color='black', transform=ccrs.PlateCarree())
    axs[2].scatter(pos5.lon, pos5.lat, color='black', transform=ccrs.PlateCarree())
    #extent = axs[2].get_extent()
    axs[0].set_title('Uncertain release point\nknown velocity field')
    axs[1].set_title('Uncertain velocity field\nknown release point')
    axs[2].set_title('Uncertain release point\nand velocity field')
    
        
    label_style = {'size': 12, 'color': 'gray'}
    for i in range(3):
        axs[i].set_extent([15.5,16.7,69.9,70.8])
        #axs[i].set_extent(extent, crs=ccrs.NorthPolarStereo(central_longitude=-24.05))
        axs[i].gridlines(draw_labels=False, x_inline=False, y_inline=False, color='gray', alpha=0.5)
        axs[i].add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
        gl = axs[i].gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
        gl.xlines = gl.ylines = False
        gl.top_labels=False
        gl.right_labels=False
        gl.rotate_labels = False
        gl.xlabel_style = gl.ylabel_style = label_style
        if i in [1,2]:
            gl.left_labels=False
        #if j in [1,2,3]:
        #    gl.left_labels=False
        #axs[i].set_aspect('equal', share=True)
    plt.tight_layout()
    plt.savefig('figures/particle_uncertainty.pdf')


plot_all()
    
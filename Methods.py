import os
import xarray as xr
import numpy as np
import re
import pandas as pd
import scipy as sp
import scipy.ndimage as ndimage
import scipy.ndimage as nd

def _normalize_DataArray_no_members(ds, alcsmin, alcsmax):
    time_index = np.arange(0, len(ds.time.values))
    
    def _select_data(time):
        return ds.isel(time=time)
    
    def _normalize_selected(dds):
        return normalize_LCS(dds, alcsmin, alcsmax)
    
    dss = xr.concat((_normalize_selected(_select_data(t)) for t in time_index), dim='time')
    return dss

def _normalize_DataArray(ds, alcsmin, alcsmax):
    members = np.arange(0,len(ds.member.values))
    time_index = np.arange(0, len(ds.time.values))

    def _select_data(member, time):
        return ds.isel(member=member).isel(time=time)
    
    def _normalize_selected(dds):
        return normalize_LCS(dds, alcsmin, alcsmax)
    
    def _concat_over_time(member):
        return xr.concat((_normalize_selected(_select_data(member, t)) for t in time_index), dim='time')
    
    dss = xr.concat((_concat_over_time(m) for m in members), dim='member')
    return dss

def normalize_LCS(LCS, alcsmin, alcsmax):
    norm_lcs = (LCS.ALCS - alcsmin)/(alcsmax-alcsmin)
    LCS['ALCS'] = (('x', 'y'), np.array(norm_lcs))
    return LCS


def mask_no_file(LCS, perc_low = 7.3, perc_high = 76.5, perc_low_custom=0, perc_high_custom=100, iter=10):
    ALCS = np.array(LCS['ALCS'])
    ALCS[np.isnan(ALCS)] = 0
    highmask = ALCS > np.percentile(ALCS, perc_high)
    lowmask = ALCS < np.percentile(ALCS, perc_low)
    ALCS[highmask] = 0
    ALCS[lowmask] = 0
    landmask = ndimage.binary_opening(ALCS, iterations=iter)
    ALCS[landmask==False] = 0
    lowmask = ALCS < np.percentile(ALCS, perc_low_custom)
    highmask = ALCS > np.percentile(ALCS, perc_high_custom)
    ALCS[lowmask] = 0
    ALCS[highmask] = 0
    
    ALCS[ALCS==0] = np.nan

    LCS['ALCS'] = (('x', 'y'), ALCS)
    
    return LCS

def mask(LCS, file, perc_low = 7.3, perc_high = 76.5, perc_low_custom=0, perc_high_custom=100, iter=10, dg = False):
    ds = xr.open_dataset(file)
    
    x0, y0 = np.array(ds.lon), np.array(ds.lat)
    nx, ny = x0.shape[0], x0.shape[1]
    x1 = np.reshape(np.array(ds.nlon), (nx, ny))
    y1 = np.reshape(np.array(ds.nlat), (nx, ny))
    ds = np.sqrt((x1-x0)**2 + (y1-y0)**2)
    ds = ds[::-1, ::-1]
    ds[np.isnan(ds)] = 0
    landmask = ndimage.binary_dilation(ds == 0)
    ALCS = np.array(LCS['ALCS'])
    ALCS[landmask] = 0
    highmask = ALCS > np.percentile(ALCS, perc_high)
    lowmask = ALCS < np.percentile(ALCS, perc_low)
    ALCS[highmask] = 0
    ALCS[lowmask] = 0
    if dg is False:
        landmask = ndimage.binary_opening(ALCS, iterations=iter)
        ALCS[landmask==False] = 0
    lowmask = ALCS < np.percentile(ALCS, perc_low_custom)
    highmask = ALCS > np.percentile(ALCS, perc_high_custom)
    ALCS[lowmask] = 0
    ALCS[highmask] = 0

    ALCS[ALCS==0] = np.nan

    LCS['ALCS'] = (('x', 'y'), ALCS)

    return LCS

def _set_up_DataArray(members=24, hours=24):
    path = os.getcwd()
    path_to_files = os.path.abspath(os.path.join(path, os.pardir))+'/scripts/LCS_files'
    mem_dirs = os.listdir(path_to_files)
    #mem_dirs.remove('DoubleGyre')
    mem_dirs.sort(key=lambda f: int(re.sub('\D*', '', f)))

    #sorting the files
    files = [f for f in os.listdir(f'{path_to_files}/{mem_dirs[0]}') if re.search('.*\_LCS\.nc', f)]
    times = [re.sub('\_h', '-', f[0:-7]) for f in files]
    times = [f'{t[0:-2]}{t[-2:]}00' if len(re.sub('(?<=^).*\-', '', t)) == 2 
                else f'{t[0:-1]}0{t[-1:]}00' for t in times]
    times, files = (list(t) for t in zip(*sorted(zip(times, files))))

    def _open_dataset(member, file):
        d = xr.open_dataset(f'{path_to_files}/{mem_dirs[member]}/{file}')
        d = d.assign(lcs_file=f'{path_to_files}/{mem_dirs[member]}/{file}', advected_file=f'{path_to_files}/{mem_dirs[member]}/{file[:-7]}.nc')
        return d

    def _concat_dataset_time(member):
        return xr.concat((_open_dataset(member, file) for file in files[0:int(hours)]), pd.Index(times[0:int(hours)], name='time'))

    ds = xr.concat((_concat_dataset_time(member) for member in range(members)), dim='member')
    return ds

def _set_up_DataArray2(members=24, hours=24, lcspath='/scripts/logs/ALCS_files'):
    '''
        Creates a DataArray using the files in LCS_files_fixed.
        The DataArray is structured by members and times as their dimensions. 
        The name of the LCS file and the advected file is also added to the DataArray.
    Args:
        members     [int]                   :   number of members to add, default = 24 (all members in Barents-2.5 EPS)
        hours       [int] or [np.ndarray]   :   The times to add to the DataArray. If [int], all times from 0 to hours will be added. If [np.ndarray] only the times in the array will be added. 
    Returns:
        ds          [DataArray]             :   an xarray DataArray containing the LCS dataset
        
    '''
    #path = os.getcwd()
    #path_to_files = os.path.abspath(os.path.join(path, os.pardir))+lcspath
    path_to_files = lcspath
    mem_dirs = os.listdir(path_to_files)
    #mem_dirs.remove('DoubleGyre')
    mem_dirs.sort(key=lambda f: int(re.sub('\D*', '', f)))

    #sorting the files
    files = [f for f in os.listdir(f'{path_to_files}/{mem_dirs[0]}') if re.search('.*\_LCS\.nc', f)]
    times = [re.sub('\_h', '-', f[0:-7]) for f in files]
    times = [f'{t[0:-2]}{t[-2:]}00' if len(re.sub('(?<=^).*\-', '', t)) == 2 
                else f'{t[0:-1]}0{t[-1:]}00' for t in times]
    times, files = (list(t) for t in zip(*sorted(zip(times, files))))

    def _open_dataset(member, file):
        d = xr.open_dataset(f'{path_to_files}/{mem_dirs[member]}/{file}')
        d = d.assign(lcs_file=f'{path_to_files}/{mem_dirs[member]}/{file}', advected_file=f'{path_to_files}/{mem_dirs[member]}/{file[:-7]}.nc')
        return d

    def _concat_dataset_time_integer(member):
        return xr.concat((_open_dataset(member, file) for file in files[0:int(hours)]), pd.Index(times[0:int(hours)], name='time'))
    
    def _concat_dataset_time_array(member):
        return xr.concat((_open_dataset(member, file) for file in _files), pd.Index(_times, name='time'))
    
    if isinstance(hours, int):
        print('instance: int')
        ds = xr.concat((_concat_dataset_time_integer(member) for member in range(members)), dim='member')
    elif isinstance(hours, np.ndarray):
        print('instance: array')
        _files = [files[index] for index in hours]
        _times = [times[index] for index in hours]
        ds = xr.concat((_concat_dataset_time_array(member) for member in range(members)), dim='member')
    return ds

def dataset_converter(x):
    x['ALCS'] = (('x', 'y'), np.log(np.sqrt(np.e**(x['ALCS'].values)))/2)
    return x

def _set_up_DataArray_masked(ds, perc_low = 7.3, perc_high = 76.5, perc_low_custom=0, perc_high_custom=100, iter=10, dg = False):
    #take in the ds created by _set_up_DataArray() and mask all ALCS fields contained in it
    members = np.arange(0,len(ds.member.values))
    time_index = np.arange(0,len(ds.time.values))

    def _select_data(member, time):
        return ds.isel(member=member).isel(time=time)
    
    def _mask_selected(dds):
        #return mask_no_file(dds, os.path.abspath(str(dds.advected_file.values)), perc_low=perc_low, perc_high=perc_high, perc_high_custom=perc_high_custom, perc_low_custom=perc_low_custom)
        return mask_no_file(dds, perc_low=perc_low, perc_high=perc_high, perc_high_custom=perc_high_custom, perc_low_custom=perc_low_custom)
    
    def _concat_over_time(member):
        return xr.concat((_mask_selected(_select_data(member, t)) for t in time_index), dim='time')

    dss = xr.concat((_concat_over_time(m) for m in members), dim='member')
    return dss

def find_correlations(alcs, topo, vel, output=True):
    topo = np.array(topo.values)
    vel = np.array(vel.values)
    x = np.arange(0,len(topo),1)
    xn = np.arange(0, len(topo), 1/2.5)
    topo = np.interp(xn, x, topo)
    vel = np.interp(xn, x, vel)

    corrt = np.corrcoef(np.array(alcs.ALCS.values)[1:], topo[1:])
    corrv = np.corrcoef(np.array(alcs.ALCS.values)[1:], vel[1:])
    corrtv = np.corrcoef(vel, topo)
    
    if output:
        print('Correlation with topography: ', corrt[0,1])
        print('Correlation with velocity: ', corrv[0,1])
        print('Correlation between topography and velocity: ', corrtv[0,1])
    return corrt[0,1], corrv[0,1], corrtv[0,1]

def stash_for_plots(axs, d):
    import matplotlib as plt
    import cartopy
    import cartopy.crs as ccrs
    
    axs.contour(d.lon_rho, d.lat_rho, d.h, transform=ccrs.PlateCarree(), colors='black', levels=10, linestyles='dashed', linewidths = 0.7)
    axs.gridlines(draw_labels=False)
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')


def DGparticles():
    from opendrift.models.oceandrift import OceanDrift
    from opendrift.readers import reader_double_gyre
    o = OceanDrift(loglevel=30)  # Set loglevel to 0 for debug information
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('drift:advection_scheme', 'runge-kutta4')

    double_gyre = reader_double_gyre.Reader(epsilon=0.25, omega=0.682, A=0.1)

    o.add_reader(double_gyre)

    x = [.85]
    y = [.9]
    lon, lat = double_gyre.xy2lonlat(x, y)
    time_step = 0.05
    from datetime import timedelta
    o.seed_elements(lon, lat, radius=.1, number=5000,
                    time=double_gyre.initial_time+timedelta(seconds=0))
    o.run(duration=timedelta(seconds=20), time_step=time_step, outfile='particle_advection7.nc')

def EPSParticles():
    from datetime import timedelta
    from opendrift.readers import reader_netCDF_CF_generic
    from opendrift.models.oceandrift import OceanDrift
    o = OceanDrift(loglevel=20)
    o.set_config('drift:advection_scheme', 'runge-kutta4')
    r = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps/barents_eps_20221001T00Z.nc', ensemble_member=0)
    o.add_reader(r)
    o.seed_elements(lon=11.4, lat=67.1, number=500000, radius=10000,
                time=r.start_time)
    o.run(duration=timedelta(hours=2), time_step=timedelta(hours=1), outfile='eps_particles3.nc')

def rotate_vectorfield(U,V,alpha):
        '''rotate wind vectors clockwise. alpha may be a scalar or an array
        alpha is in degrees
        returns u,v '''
        alpha = np.array(alpha)*sp.pi/180
        alpha = alpha.flatten()
        R = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)] ])
        shpe = U.shape
        origwind = np.array((U.flatten(), V.flatten()))
        if len(R.shape)==2:
            rotwind = dot(R, origwind) # for constant rotation angle
        else:
            # for rotation angle given as array with same dimensions as U and V:
            # k-loop with rotwind(k) = dot(R(i,j,k), origwind(j,k)) (einstein summation indices)
            rotwind = np.einsum("ijk,ik -> jk", R, origwind)  # einstein summation indices
        Urot, Vrot = rotwind[0,:], rotwind[1,:]
        Urot = Urot.reshape(shpe)
        Vrot = Vrot.reshape(shpe)
        return Urot, Vrot
    
def north_direction(lat):
    '''get the north direction relative to image positive y coordinate'''
    dlatdx = nd.filters.sobel(lat,axis=1,mode='constant',cval=sp.nan) #gradient in x-direction
    dlatdy = nd.filters.sobel(lat,axis=0,mode='constant',cval=sp.nan)
    ydir = lat[-1,0] -lat[0,0] # check if latitude is ascending or descending in y axis
    # same step might have to be done with x direction.
    return np.arctan2(dlatdx,dlatdy*np.sign(ydir) )*180/sp.pi

if __name__ == '__main__':
    DGparticles()


    



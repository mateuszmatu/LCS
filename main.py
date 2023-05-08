from ParticleAdvector import BackwardsAdvection
from LCS import FTLE
import xarray as xr
import numpy as np
from joblib import Parallel, delayed

def _mkdir(member):
    import os
    if not os.path.exists(store_file):
        os.mkdir(store_file)
    if not os.path.exists(f'{store_file}/member{member}'):
        os.mkdir(f'{store_file}/member{member}')

def CreateLCSField(lons, lats, ts, sep, dur, date, member, output, at_time=0):
    _mkdir(member)
    advector = BackwardsAdvection(lons, lats, ts, sep, dur, date, at_time)
    outfile = advector.displace_one_member(member, output)
    LCS = xr.open_dataset(FTLE(f'{output}.nc', f'{output}_LCS'))

def Run(i, date, at_time=24):
    output=f'{store_file}/member{i}/{date}_h{at_time}'
    CreateLCSField(lons=[7,23], lats=[67,69.9], ts=-900, sep=1000, dur=2, date=date, member=i, output=output, at_time=at_time)

if __name__ == '__main__':
    """
    import os
    for i in range(24):
        os.mkdir(f'logs/ALCS_files_2h_april/member{i}')
    """
    store_file = 'ALCS_files_2h_april'
    import sys
    sel = int(sys.argv[1])
    times=np.arange(0,24,2)

    dates = [f'202204{d:02d}' for d in range(1,5)]
    #dates = ['20221002', '20221003']
    curr_date = dates[sel-1]
    for t in times:
        for i in range(24):
            Run(i, curr_date, at_time=int(t))
    
    
    """
    for t in times:
        for i in range(24):
            Run(i, int(t), curr_date)
    """
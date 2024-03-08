import os
import re
from datetime import datetime
from datetime import timedelta

def files_from_thredds(date, file):
    '''
        Find the threeds files from {file}.txt from the previous day of the specified date.
    Args:
        date    [str]   :   the date from which files should be found
        file    [str]   :   file with the urls
    Returns:
        arr     [list]  :   a list of the file names from the specified date
    '''
    odate = datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]))
    ndate = odate-timedelta(days=1)
    date = re.sub('-', '', str(ndate)[0:10])
    date_regex = rf'{date}'
    arr = []
    with open(file) as file:
        lines = file.readlines()
        for line in lines:
            if re.findall(date_regex, line):
                line = re.sub('\n', '', line)
                arr.append(line)
    arr.sort()
    if len(arr) == 4:
        pass
    else:
        raise ValueError('Something happened, list of files does not contain all members.')

    return arr

def files_from_lustre(date):
    path = os.path.abspath('/lustre/storeB/project/fou/hi/barents_eps/eps')
    files = os.listdir(path)
    arr = []
    for file in files:
        if re.findall(date, file) and re.findall('eps', file):
            arr.append(path+'/'+file)
    arr.sort()

    if len(arr) == 4:
        pass
    else:
        raise ValueError('Something happened, list of files does not contain all members.')

    return arr

def correct_file(arr, member):
    '''
        Finds the correct file in the list of files from files_from_thredds(), based on which member is chosen
    Args:
        arr     [list]  :   a list of four filenames
        member  [int]   :   the member to look at
    Returns:
        file    [str]   :   a string of the correct filename
        _member [int]   :   an integer corresponding to the requested member, but in the selected file
    '''
    if member >= 0 and member < 6:
        file = arr[0]
        _member = member
    elif member >= 6 and member < 12:
        file = arr[1]
        _member = member - 6
    elif member >= 12 and member < 18:
        file = arr[2]
        _member = member - 12
    elif member >= 18 and member < 24:
        file = arr[3]
        _member = member - 18
    return file, _member

def hour_adjustment(member):
    '''
        Adds hours to the file based one which member it is
    Args: 
        member  [int]   :   the member used
    Returns:
        addhour [int]   :   hours that are supposed to be added to the file
    '''
    if member >= 0 and member < 6:
        addhour = 24
    elif member >= 6 and member < 12:
        addhour = 18
    elif member >= 12 and member < 18:
        addhour = 12
    elif member >= 18 and member < 24:
        addhour = 6
    return addhour

def name_from_lon_lat(lon,lat):
    '''
        Creates a string name created from given longitude and latitude
    Args:
        lon     [list]  :   a lsist of min and max lon
        lat     [list]  :   a list of min and max lat
    Returns:
        name    [str]   :   a string with lons and lats
    '''
    name = f'lon{str(lon[0])}_{str(lon[1])}_lat{str(lat[0])}_{str(lat[1])}'
    name = re.sub('\.','-',name)
    return name

if __name__ == '__main__':
    arr = files_from_lustre('20230731')



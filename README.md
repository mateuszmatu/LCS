![](https://github.com/mateuszmatu/LCS/blob/master/gifs/DG_animation.gif)
# Todo
- [ ] Test on more data to check for robustness
- [ ] Implement 3D LCSs
- [ ] Implement a method for detecting "true" LCSs, following Farazmand and Haller 2012

# Important: This is still work in progress
# Any feedback or suggestions are appreciated

This is a repository containing more general scripts for computing hyperbolic LCSs, utilizing the FTLE approach.
The goal of the files in this folder is to make the scripts as robust as possible, so that they can be used with any model in any region. 
These scripts compute LCSs in the ocean, but should work for 2D atmosphere as well (not tested).

# Dependencies

Particles are advected using MET Norways OpenDrift python module. This module and an installation guide can be found on https://opendrift.github.io/.
Other packages needed can all be installed normally through 
```
pip install [package name]
```

# Using the software

The goal of this software is for it to be simple to use, and yet general. Particles are initiated by the ```Advection``` class found in ```ParticleAdvector.py```, and advected with ```Advection.run()```. This returns a DataArray containing the initial and final particle positions.  
Example:
```
t = Advection(file, lons, lats, ts, sep, dur, start)
parts = t.run()
```
[file] is the name of the file contaning the velocity field, [lons] is a list containing min and max longitudes, [lats] is a list containing min and max latitudes, [ts] is the time step for the integration in seconds, [sep] is the initial separation between particles in meters, [dur] is the simulation duration in hours and [start] is a datetime object containing the simulation start time. 
If the velocity fields come from an ensemble model, the specific ensemble member can be specified with ```t.run(ensemble_member=....)```.

Once the particles have been advected, they can be provided to the ```FTLE``` function in ```LCS.py```.
```
LCS = FTLE(parts, example_model_file)
```
Where [example_model_file] is again the name of the file containing the velocity fields. It is important to provide this, as the software conducts a coordinate transformation based on the datasets projection. 

# Final note
The sign of [ts] determines whether attracting (negative ts) or repelling (positive ts) LCSs are computed. Note that the software by default assumes that attracting LCSs are computed. For repelling LCSs, the outputted LCSs have to be flipped
```
LCS = LCS[::-1,::-1]
```
The reason for this is still unknown.

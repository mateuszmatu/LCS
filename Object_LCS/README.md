![](https://github.com/mateuszmatu/LCS/blob/master/gifs/DG_animation.gif)

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
parts = t.run()#(ensemble_member=0)
```


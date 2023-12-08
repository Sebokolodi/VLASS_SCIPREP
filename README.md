# VLASS_SCIPREP
Scripts to use for VLASS (polarization) science pre-processing pipeline. The available scripts are 
  1. conv_align.py: This does two things: i) it re-centers the crval1/2 of the cubes which got misaligned
     during image-plane correction for w-terms (see https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide)
     for details. Secondly, it convolves the cubes to user specified resolution.
  2. stokes_freq_cube.py: This script merges the individual spectral windows into a frequency cube, and splits the
     the cubes based on Stokes parameters I, Q and U. So we went up with 3 sets of data cubes (each per Stokes)
     with frequency as a third axis.
  3. vlass_cutout_3D.py: Creates cutouts for both the mfs and cubes.


## Dependancies:
1. CASA: needed for conv_align and stokes_freq_cube. See installation options here: https://casadocs.readthedocs.io/en/v6.2.0/notebooks/usingcasa.html.
2. Montage: needed for vlass_cutout_3D.
3. Astropy: needed for stokes_freq_cube and vlass_cutout_3D scripts.
4. Pandas: needed for vlass_cutout_3D.
5. Reproject: needed for vlass_cutout_3D.
   
Note: Running astropy inside CASA is not as straightforward esp when using full installation version of CASA.
You can install astropy inside CASA by following: https://stackoverflow.com/questions/52289107/installation-of-astropy-in-casa (install pip) and
https://astropy-cjhang.readthedocs.io/en/latest/install.html (install astropy).

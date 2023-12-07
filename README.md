# VLASS_SCIPREP
Scripts to use for VLASS (polarization) science pre-processing pipeline. The available scripts are 
  1. conv_align.py: This does two things: i) it re-centers the crval1/2 of the cubes which got misaligned
     during image-plane correction for w-terms (see https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide)
     for details. Secondly, it convolves the cubes to user specified resolution.
  2. stokes_freq_cube.py: This script merges the individual spectral windows into a frequency cube, and splits the
     the cubes based on Stokes parameters I, Q and U. So we went up with 3 sets of data cubes (each per Stokes)
     with frequency as a third axis.
  3. vlass_cutout_3D.py: Creates cutouts for both the mfs and cubes. 

# VLASS_SCIPREP
Scripts to use for VLASS (polarization) science pre-processing pipeline. The available scripts are 
  1. conv_align.py: This does two things: i) it re-centers the crval1/2 of the cubes which got misaligned
     during image-plane correction for w-terms (see https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide)
     for details. Secondly, it convolves the cubes to user specified resolution.
  2. stokes_freq_cube.py: This script merges the individual spectral windows into a frequency cube, and splits the
     the cubes based on Stokes parameters I, Q and U. So we went up with 3 sets of data cubes (each per Stokes)
     with frequency as a third axis.
  3. vlass_cutout_3D.py: Creates cutouts for both the mfs and cubes.


## Dependancies
1. CASA: needed for conv_align and stokes_freq_cube. See installation options here: https://casadocs.readthedocs.io/en/v6.2.0/notebooks/usingcasa.html.
2. Montage: needed for vlass_cutout_3D.
3. Astropy: needed for stokes_freq_cube and vlass_cutout_3D scripts.
4. Pandas: needed for vlass_cutout_3D.
5. Reproject: needed for vlass_cutout_3D.
   
Note: Running astropy inside CASA is not as straightforward esp when using full installation version of CASA.
You can install astropy inside CASA by following: https://stackoverflow.com/questions/52289107/installation-of-astropy-in-casa (install pip) and
https://astropy-cjhang.readthedocs.io/en/latest/install.html (install astropy).

## Cutout Script
The cutout script takes as input the source ra and dec [in degrees] of interest, the epoch e.g. 1, 2, or 3, and cutout size in arcmins. You need
to also specify a path to directory with fits files to cut into small chuncks, and path to directory to store both the cutout images and mosaic image (if requested).

Though not user to specify, the script requires Subtile information. This is store inside the folder subtile_catalogs. Each epoch (including both halves)
has its own subtile information and is obtained from the CIRADA continuum group. The idea is that we will store all the subtiles csv files inside
this folder and access them when we call the script. At the moment (7th of December), we only have info about epoch 1 and 2.  The stored 
csv filename should take the format 
    
     CIRADA_VLASS2Q*_subtile.csv
    
where 2 can be replaced by 1, 2, or 3. So far we do not have epoch 3. 

This subtile csv file contain important parameters such as crval1, crval2, latpole,
lonpole, bmaj, bmin, and bpa used for generating the images. We use this information to generate template headers for each subtile. Other information that 
we use that aren't available inside the subtile csv file include: the naxis and cdelt1/2. These are confirmed to be the same for all of the VLASS image data. 
So are fixed into the code. With the template header, we are able to determine the boundaries of the images, and thus, determine which subtile a source of interest 
appears. Since the VLASS images overlap to some extend at the edges, a source can appear in more than one subtile.



## Pending: To dos

1. Write a sript that replaces missing spectral windows with channels containing nans. This should be incorporated before converting to frequency cubes,
   so that the header information is correct (at this point, jumps in frequency due to missing SPW results in inaccurate frequency range).
2. To include a portion inside vlass_cutout_3D that considers overlaps when the requested cutout size extents into another subtile. In which case, the user
   has the option to request mosaicking.
3.  





# VLASS_SCIPREP
Scripts to use for VLASS (polarization) science pre-processing pipeline. The available scripts are 
  1. conv_align.py: Does two things: i) it re-centers the crval1/2 to the appropriate values (see https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide
     for details). Secondly, it convolves the cubes to user specified resolution.
  2. stokes_freq_cube.py: Merges the individual spectral windows into frequency cubes, and splits the
     the cubes into separate Stokes parameters, namely Stokes I, Q and U. 
  3. vlass_cutout_3D.py: Creates cutouts for both the mfs and cubes.


## Dependancies
1. CASA: Needed for conv_align and stokes_freq_cube. See installation options here: https://casadocs.readthedocs.io/en/v6.2.0/notebooks/usingcasa.html.
2. Montage: Needed for vlass_cutout_3D.
3. Astropy: Needed for stokes_freq_cube and vlass_cutout_3D scripts.
4. Pandas: Needed for vlass_cutout_3D.
5. Reproject: Needed for vlass_cutout_3D.
   
Note: Running astropy inside CASA is not as straightforward esp when using full installation version of CASA.
You can install astropy inside CASA by following: https://stackoverflow.com/questions/52289107/installation-of-astropy-in-casa (install pip) and
https://astropy-cjhang.readthedocs.io/en/latest/install.html (install astropy).

## Cutout Script
The cutout script requires the source RA and DEC (in degrees), the epoch (e.g., 1, 2, or 3), and the cutout size (in arcminutes) as input. 
Additionally, you need to provide the path to the directory containing the FITS files (MFS or frequency cubes) and the path to the directory
where the cutout images and mosaic image (if requested) will be stored. Finally, you can specify whether you want to create a mosaic.

Although the user does not need to specify it, the script also takes in subtile information stored in the "subtile_catalogs" folder. 
Each epoch has its own subtile information (e.g. 1.1 and 1.2), obtained from the CIRADA continuum group. The idea is to store all the subtiles CSV files
inside this folder and access them when the script is called. The stored CSV filename should following format:
    
     CIRADA_VLASS2Q*_subtile.csv
    
where 2 can be replaced by 1, 2, or 3. So far we do not have epoch 3. 

The subtile CSV file contains crucial parameters like crval1, crval2, latpole, lonpole, bmaj, bmin, and bpa from the original images. 
This information allows us to create template headers for each subtile, eliminating the need to download images to retrieve this info. 
Additional info not in the subtile CSV include naxis and cdelt1/2, are consistent across all VLASS image data and are fixed in the code. 
With the template header, we can determine the image boundaries and identify sources belonging to a subtile. Also, the VLASS images overlap at the edges, 
as such, some sources appear in more than one subtile. For the latter, there is a mosaic option to put together subtiles overlapping at a source location 
(yet to be implemented). 



## Pending: To dos

1. Write a sript that replaces missing spectral windows with channels containing nans. This should be incorporated before converting to frequency cubes,
   so that the header information is correct (at this point, jumps in frequency due to missing SPW results in inaccurate frequency range).
2. Option to mosaic overlapping subtiles.
   





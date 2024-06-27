#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# VLASS 3D cutout strategy (we can extend it to POSSUM later).

# User inputs: 1) source coordinates (ra & dec).
#              2) Epoch (e.g. 1, 2, 3).
#              3) Cutout size in arcmin.
#              4) Mosaic or not. 

# Data:    1) MFS and cubes from the folder. (Ideally we want to acccess these data online/archive).
#          2) CIRADA subtile catalog: which contains subtile name, crval, latpole/lonpole. 

# Processing: 1) Take a coordinate and Epoch to search for fields to which these appear.  
#             2) Return the field over a specific cutout.
#             3) If no mosaic, return all the cutouts for which the source appears. 
#                If mosaic specified, then mosaic all the fields to which the source appears.            

# Validation: 1) Ensure that the peak source is the same for mosaicked and non-mosaicked.


# Test: 1) on mfs image.
#       2) on cube images.


import glob
import pandas as pd
import numpy
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
import montage_wrapper as montage
from MontagePy.main import mAddCube, mViewer
import tempfile
import os, sys 
import shutil, io
from reproject import reproject_interp
import argparse



def reference_header(naxis, cdelt, crval1, crval2, latpole, lonpole, bmaj, bmin, bpa):
    
    """ Create a reference header used to determine overlaping subtiles. """

    cdelt = abs(cdelt) 
    
    hdr  = "SIMPLE  =                    T / file does conform to FITS standard \n"
    hdr += "BITPIX  =                  -32 / number of bits per data pixel \n"
    hdr += "NAXIS   =                    2 / number of data axes  \n"
    hdr += "NAXIS1  =                %d / length of data axis 1 \n" %naxis
    hdr += "NAXIS2  =                %d / length of data axis 2 \n" %naxis
    hdr += "EXTEND  =                    T / No FITS extensions are present \n"
    hdr += "BMAJ    =                %r \n"%bmaj
    hdr += "BMIN    =                %r \n"%bmin
    hdr += "BPA     =                %r \n"%bpa
    hdr += "LONPOLE =                %r \n"%lonpole                                                  
    hdr += "LATPOLE =                %r \n"%latpole                                                
    hdr += "PC1_1   =   1.000000000000E+00 \n"                                                  
    hdr += "PC2_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC3_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC4_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC1_2   =   0.000000000000E+00 \n"                                                  
    hdr += "PC2_2   =   1.000000000000E+00 \n"                                                  
    hdr += "PC3_2   =   0.000000000000E+00 \n"                                                  
    hdr += "PC4_2   =   0.000000000000E+00 \n"                                                 
    hdr += "PC1_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC2_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC3_3   =   1.000000000000E+00 \n"                                                 
    hdr += "PC4_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC1_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC2_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC3_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC4_4   =   1.000000000000E+00 \n"                                                 
    hdr += "CTYPE1  = 'RA---SIN'   \n"                                                         
    hdr += "CRVAL1  =   %r \n"%crval1                                                
    hdr += "CDELT1  =  %r  \n"%(-cdelt)                                                 
    hdr += "CRPIX1  =  %r \n"%(naxis/2.0)                                                  
    hdr += "CUNIT1  = 'deg      '  \n"                                                       
    hdr += "CTYPE2  = 'DEC--SIN' \n"                                                            
    hdr += "CRVAL2  =  %r \n"%crval2                                                  
    hdr += "CDELT2  =  %r \n"%cdelt                                                 
    hdr += "CRPIX2  =  %r \n"%(naxis/2.)                                                 
    hdr += "CUNIT2  = 'deg     ' \n"  
                    
    return hdr


def read_csv_subtiles(csv_file):
    
    """ Read csv file"""
    
    data = pd.read_csv(csv_file)
    
    return data
    

def find_nearest_subtiles(source_ra, source_dec, subtile_crval1, subtile_crval2):
    
    """From the subtile csv, determine the subtiles nearest to the source coordinates
    
    source_ra: input source ra in degrees.
    source_dec: input source dec in degrees.
    
    subtile_crval1/2: crval values of the subtiles (as read from the subtile 
                      csv catalog.  We get this catalog from the CIRADA continuum).
                      
    Returns
    
    Index corresponding to the 6 nearest subtiles. 
    """
    
    ref_pos  = SkyCoord(ra = source_ra*u.degree, dec=source_dec*u.degree, frame='fk5')
    subtile_pos = SkyCoord(ra=subtile_crval1*u.degree, dec=subtile_crval2*u.degree, frame='fk5')  
    sep = ref_pos.separation(subtile_pos).deg

    nearest_neighbours = sorted(sep)[:6]  # take the 6 nearest subtile.
    nearby_index = [] # return the subtile index.
    for near in nearest_neighbours:
        nearby_index.append(numpy.where(sep == near)[0][0])

    return nearby_index
    


def in_or_out_field(ra, dec, wcs):
    
    # convert ra and dec in deg to pixel coordinate.
    rapix, decpix = wcs.wcs_world2pix(ra, dec, 0)

    return rapix, decpix

    
def overlap_subtile(source_ra, source_dec, nearby_subtile_data):

    """ Determines the subtile(s) for which the source appears. """

    # Considering the nearest subtiles.   
    index = []    

    tile = [d for d in nearby_subtile_data['Tile']]
    subtile = [d for d in nearby_subtile_data['Subtile']]
    bmaj = [d for d in nearby_subtile_data['BMAJ']]
    bmin = [d for d in nearby_subtile_data['BMIN']] 
    bpa  = [d for d in nearby_subtile_data['BPA']]
    latpole = [d for d in nearby_subtile_data['LATPOLE']]
    lonpole = [d for d in nearby_subtile_data['LONPOLE']]
    crval1 = [d for d in nearby_subtile_data['CRVAL1']]
    crval2 = [d for d in nearby_subtile_data['CRVAL2']]
    
    for i in range(len(subtile)):
        hdr = reference_header(naxis=6202, cdelt=1.666666666667E-04, crval1=crval1[i],
                          crval2=crval2[i], latpole=latpole[i], lonpole=lonpole[i], 
                          bmaj=bmaj[i], bmin=bmin[i], bpa=bpa[i])
        hdr = fits.Header.fromstring("""%s"""%hdr, sep='\n')
        wcs = WCS(hdr, naxis=2)
        
        rapix, decpix = in_or_out_field(source_ra, source_dec, wcs)
        naxis1 = hdr['naxis1']
        naxis2 = hdr['naxis2']
  
        # If the pix coordinate is within naxis and is positive, then the coordinate
        # exists inside the image. 
        # For a cutout, we will need to use rapix = rapix -/+ cutout_size/2, same for dec too.
        # To determine if the cutout boundary is within an image. 
        
        if (rapix < naxis1) & (decpix < naxis2) & (rapix > 0) & (decpix > 0):
            index.append(i)  
               
    return nearby_subtile_data.iloc[index]
    

def read_image(image):
    
    """Reads fits data using astropy, and returns a 2D plane data,
    header and wcs. If you provide a cube, it will use the first plane."""
                                                                              
    with fits.open(image) as hdu:
        #hdr = hdu[0].header
        data = hdu[0].data
        #wcs = WCS(hdr).celestial
    imslice = numpy.zeros(data.ndim, dtype=int).tolist()
    imslice[-1] = slice(None)
    imslice[-2] = slice(None)
    return data[tuple(imslice)]#, hdr, wcs
    
  
def cutout_ref_header(naxis, crval1, crval2, img_header, ref_header):
    
    """ Create a reference header for reprojecting cutouts to the same 
        grid and reference coordinates. 
        
    naxis: the size of the cutout. Value used for both naxis 1 and 2. 
    crval1/2: are the source ra and dec. These form the central coordinates
              of a cutout.
    img_header: Reads header information of the subtile in question.
    ref_header: A reference header used to determine the reference lonpole, latpole. 
                naxis3 (freq), and naxis 4 (stokes).
                
    Returns
    
    Header to use for reprojecting.
      
    """
    cdelt = abs(ref_header['cdelt1'])
    
    
    hdr  = "SIMPLE  =                    T / file does conform to FITS standard \n"
    hdr += "BITPIX  =                  -32 / number of bits per data pixel \n"
    hdr += "NAXIS   =                    4 / number of data axes  \n"
    hdr += "NAXIS1  =                %d / length of data axis 1 \n" %naxis
    hdr += "NAXIS2  =                %d / length of data axis 2 \n" %naxis
    hdr += "NAXIS3  =                %d / length of data axis 3 \n" %ref_header['naxis3']
    hdr += "NAXIS4  =                %d / length of data axis 4 \n" %ref_header['naxis4']
    hdr += "EXTEND  =                    T / No FITS extensions are present \n"
    hdr += "LONPOLE =   %r \n"%ref_header['lonpole']                                                  
    hdr += "LATPOLE =   %r \n"%ref_header['latpole'] 
    hdr += "BMAJ    =   %r \n"%img_header['bmaj']                                               
    hdr += "BMIN    =   %r \n"%img_header['bmin']                                                 
    hdr += "BPA     =   %r \n"%img_header['bpa'] 
    hdr += "OBJECT  =   %r \n"%img_header['object']
    hdr += "BTYPE   =   %r \n"%img_header['BTYPE']                                                                                                           
    hdr += "BUNIT   =   %r  /Brightness (pixel) unit  \n"%img_header['BUNIT']                                                                                   
    hdr += "RADESYS =   %r     \n"%img_header['RADESYS']      
    hdr += "PC1_1   =   1.000000000000E+00 \n"                                                  
    hdr += "PC2_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC3_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC4_1   =   0.000000000000E+00 \n"                                                  
    hdr += "PC1_2   =   0.000000000000E+00 \n"                                                  
    hdr += "PC2_2   =   1.000000000000E+00 \n"                                                  
    hdr += "PC3_2   =   0.000000000000E+00 \n"                                                  
    hdr += "PC4_2   =   0.000000000000E+00 \n"                                                 
    hdr += "PC1_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC2_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC3_3   =   1.000000000000E+00 \n"                                                 
    hdr += "PC4_3   =   0.000000000000E+00 \n"                                                 
    hdr += "PC1_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC2_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC3_4   =   0.000000000000E+00 \n"                                                 
    hdr += "PC4_4   =   1.000000000000E+00 \n"                                                 
    hdr += "CTYPE1  =   'RA---SIN'   \n"                                                         
    hdr += "CRVAL1  =   %r \n"%crval1                                                
    hdr += "CDELT1  =   %r  \n"%(-cdelt)                                                 
    hdr += "CRPIX1  =   %r \n"%(naxis/2.0)                                                  
    hdr += "CUNIT1  =   'deg      '  \n"                                                       
    hdr += "CTYPE2  =   'DEC--SIN' \n"                                                            
    hdr += "CRVAL2  =   %r \n"%crval2                                                  
    hdr += "CDELT2  =   %r \n"%cdelt                                                 
    hdr += "CRPIX2  =   %r \n"%(naxis/2.)                                                 
    hdr += "CUNIT2  =   'deg     ' \n" 
    hdr += "CTYPE3  =   'FREQ    ' \n"                                                          
    hdr += "CRVAL3  =   %r \n"%ref_header['crval3']                                                 
    hdr += "CDELT3  =   %r \n"%ref_header['cdelt3']                                                   
    hdr += "CRPIX3  =   %r \n"%ref_header['crpix3']                                                  
    hdr += "CUNIT3  =   'Hz      '  \n"                                                           
    hdr += "CTYPE4  =   'STOKES  '  \n"                                                           
    hdr += "CRVAL4  =   %r \n"%ref_header['crval4']                                                   
    hdr += "CDELT4  =   %r \n"%ref_header['cdelt4']                                                  
    hdr += "CRPIX4  =   %r \n"%ref_header['crpix4']                                                  
    hdr += "CUNIT4  = '       ' \n"  
    hdr += "PV2_1   =   %r \n"%img_header['PV2_1']                                                
    hdr += "PV2_2   =   %r \n"%img_header['PV2_2']                                                 
    hdr += "RESTFRQ =   %r /Rest Frequency (Hz) \n"%img_header['RESTFRQ'] 
              
    return hdr
 
    
    
def trim_tile(tile_hdu, position, cutout_size):
    
    """Trims to cutout size around a given center 
    
     tile_hdu: Is an HDU of a subitle. It consists of both 
               the data (tile_hdu[0].data) and header (tile_hdu[0].header).
     position: source ra and dec in skyCoord. 
     cutout_size: in arcmin. 
     
     returns
     
     A cutout HDU reprojected to the reference coordinate. 
     
    """
    
    if tile_hdu is None:
        return None
    
    # read the HDU.
    tile_header = tile_hdu[0].header
    tile_wcs = WCS(tile_header)
    tile_data = tile_hdu[0].data
    naxis = tile_wcs.naxis
    
    # convert cutout_size to pixels.
    pixel_size = abs(tile_header['cdelt1'])
    size = round(cutout_size/pixel_size)
    
    # We work in 2D, even for cubes (we iterate over frequency).
    # So, if an image has naxis > 2, we remove the other axis for the purpose 
    # of reprojecting. 
    while naxis > 2:
        tile_wcs = tile_wcs.dropaxis(2)
        naxis -=1
    
    # find which axis is frequency and Stokes parameter.
    if tile_header["CTYPE3"] == "STOKES":
        stokes = tile_header["CRVAL3"]  
        num_freq = tile_data.shape[0]
        
    if tile_header["CTYPE4"] == "STOKES":
        stokes = tile_header["CRVAL4"]  
        num_freq = tile_data.shape[1]
    
    #define a reference header. This is different from the input subtile
    #because now the cutout is centered about the source_ra, source_dec.
    ref_hdr = cutout_ref_header(naxis=size,  crval1=source_ra, 
                      crval2=source_dec, ref_header=tile_header,
                      img_header=tile_header)
    ref_hdr = fits.Header.fromstring("""%s"""%ref_hdr, sep='\n')
    ref_wcs = WCS(ref_hdr)
    
    # drop other axis for a reference header.
    naxis = ref_hdr['naxis']
    while naxis > 2:
        ref_wcs = ref_wcs.dropaxis(2)
        naxis -=1
    
    # if is a cube. 
    if num_freq > 1:
        # loop over the frequency planes.
        if tile_header['CTYPE3'] == 'FREQ':   
            cube = numpy.zeros((1, num_freq, size, size), dtype='float32')
            
            for plane in  range(0, num_freq):
                #reproject and create a cutout for each plane.
                array, footprint = reproject_interp(input_data = (tile_data[0, plane, ...], tile_wcs),
                       output_projection=ref_wcs, shape_out=(size, size)) 
                cube[0, plane, ...] = array
                
        if tile_header['CTYPE4'] == 'FREQ':
            # if freq axis is the axis 4, still keep freq axis in the cutout axis 3. 
            cube = numpy.zeros((1, num_freq, size, size), dtype='float32')
            for plane in range(0, nfreq):
                array, footprint = reproject_interp(input_data = (tile_data[plane, 0, ...], tile_wcs),
                       output_projection=ref_wcs, shape_out=(size, size)) 
                cube[0, plane, ...] = array
                
        tile_cutout_data = cube
    else:
        # mfs image
        img_data = numpy.squeeze(tile_data)      
        array, footprint = reproject_interp(input_data=(img_data, tile_wcs), 
             output_projection=ref_wcs, shape_out=(size, size)) 
        tile_cutout_data = array
        #ref_hdr = ref_wcs.to_header()
                
    
    cutout_header = ref_hdr
    trimmed = fits.PrimaryHDU(tile_cutout_data, header=cutout_header)
       
    return trimmed
 
   
def mosaic_images(cutouts, outdir, outprefix, debug=True):

    """Mosaics cutout images. 
    
    cutout: A list of HDUs for cutouts to mosaic.
    outdir: the output directory to store the final mosaic image. 
    outprefix: prefix to use for output mosaic, and temporary files 
               generated for Montage mosaic such as tbl, imglist, temporary header. 
               Note that for now, the prefix is alreaday defined within the code as 
               'vlass_cutout_mosaic_epoch_source_pos' . In the future, this can be left
               for a user to specify. 
    debug: This is primary for Montage. As we generate tenmporary image list, table and
           header, it best to print these so that we see when something has gone wrong. 
           
    returns:
    
    It outputs the mosaic files. 
    
    """
    
    #Create a temporary folder to store cutouts images needing mosaic.
    #In future, we want to use the actual tempfile format.
  
    temporary_indir = 'store_cutouts/' 
    if not os.path.isdir(temporary_indir):
         os.mkdir(temporary_indir)
         
    imglst_file = temporary_indir+ 'imglist_%s.list'%outprefix
    imgtab_file = temporary_indir+ 'imgtable_%s.tbl'%outprefix
    imghdr_file = temporary_indir+ 'imgheader_%s.hdr'%outprefix
    
    outimage = outdir + outprefix + '.fits'  

    with open(imglst_file, 'w') as f:
        fitsname = []
        for i, c in enumerate(cutouts):
            temp_fitsname = '{directory}/{name}.fits'.format(directory=temporary_indir,
                   name=outprefix + '_%d'%i) 
                           
            with open(temp_fitsname, 'wb') as tmp:
                img = fits.PrimaryHDU(c.data, header=c.header)
                img.writeto(tmp)
                
            fitsname.append(os.path.basename(temp_fitsname))
         
        #generate the image list for the specific images.  
        nchar = max([len(n) for n in fitsname ])
        f.writelines('| %s|\n'%('fname'.rjust(nchar)))
        f.writelines('| %s|\n'%('char'.rjust(nchar)))
        for fits_file in fitsname:
            f.write(' %s\n'%fits_file)
        f.close()
        
    #first generate the table file, from the         
    montage.mImgtbl(directory=temporary_indir, images_table=imgtab_file, 
            img_list=imglst_file, debug=debug)
     
    # to define the header, we use one of the cutout images, these are 
    # already in the proper projection, so what's left is to add them.
    # we need to test this further, to make sure that's what is happening.     
    montage.mGetHdr(in_image=temp_fitsname, img_header=imghdr_file)
    
    #get frequencies axis. This is in axis 3 as defined in the cutout cubes. 
    #An mfs might not have this axis defined. So we might get an error
    #when we try to read 'naxis3' from its header. So for now, 
    #we will try get naxis3, if it exits, and is > 1 then call maddCube
    #if is ==1 rather, call a 2D madd. 
    #if the input is an mfs without naxis3, run except and use 2D madd to
    #generate a mosaic. 
    #TODO: there might be better ways to do this. 
    try:       
        nfreq = c.header['naxis3']
        # Add the images. Store this file to the outdir specified
        if nfreq > 1: #
            mAddCube(path=temporary_indir, tblfile=imgtab_file,
                template_file=imghdr_file, outfile=outimage, coadd=0, debug=debug)
    
        else: #mosaic mfs images.
            montage.mAdd(images_table=imgtab_file, template_header=imghdr_file,
                 out_image=outimage, img_dir=temporary_indir, no_area=True)
                 
    except: #mfs image.
        montage.mAdd(images_table=imgtab_file, template_header=imghdr_file,
             out_image=outimage, img_dir=temporary_indir, no_area=True)
    
    # Then delete the temporary folder named store_cutouts.
    shutil.rmtree(temporary_indir)
    
 
def cutout_image(source_overlap_subtiles, cutout_size, fits_indir, 
        cutout_outdir, epoch=2.1, do_mosaic=False):
        
        
    """This is a main function that allows a user to create cutouts 
       and specify whether to mosaic of not.
    
    
    source_overlap_subtiles: these a table-like (csv) containing information
          about the subtiles for which the source in questions appears. 
   
    cutout_size: cut out size in arcmin. 
    
    fits_indir: a directory containg the images to cut into small chuncks.
    cutout_outdir: a directory to store the cutout images, and also the mosaic file.
    epoch: 
    do_mosaic: either true or false. If true, generate a mosaic. 
     
     
    returns
    
    cutout fits image stored inside cutout_outdir
    if mosaic is true, also mosaic cutout. 
     
    """    

    #check if the cutout_outdir exist, if not, generate one.
    if not os.path.isdir(cutout_outdir):
         print("The output directory %s is not found. Generating one..."%cutout_outdir)
         os.mkdir(cutout_outdir)

    tile = [d for d in source_overlap_subtiles['Tile']]
    subtile = [d for d in source_overlap_subtiles['Subtile']]
    print('Overlaping subtiles %s'%subtile)
        
    cutouts = []
    for i in range(len(subtile)):
    
        print(fits_indir)
        fitsimage = glob.glob(fits_indir + 'VLASS%s*%s.%s*.fits'%(epoch, tile[i], subtile[i]))
     
        if len(fitsimage) > 0:
        
            for img in fitsimage:
                print(img)
                hdu = fits.open(img)
                position = SkyCoord(source_ra, source_dec, unit='deg')
            
                if source_dec < 0:
                    source_pos = '%s%s'%(source_ra, source_dec)
                if source_dec > 0:
                    source_pos = '%s+%s'%(source_ra, source_dec)
            
                trimmed = trim_tile(tile_hdu=hdu, position=position, cutout_size=cutout_size)
                cutouts.append(trimmed) 
            
                # Define Stokes parameters. We do this in order to determine
                # whether the input image is Stokes I, Q, U and V.
               
                if trimmed.header['CTYPE4'] == "STOKES":
                    stokes = trimmed.header['CRVAL4'] 
                    print(trimmed.header['CRVAL4'])
                    if int(stokes) == 1:
                        stokesid = "i"
                    elif int(stokes) == 2:
                        stokesid = "q"
                    elif int(stokes) == 3:
                        stokesid = "u"
                    elif int(stokes) == 4:
                        stokesid = "v" 

                elif trimmed.header['CTYPE3'] == "STOKES":
                    stokes = trimmed.header['CRVAL3'] 
                    print(trimmed.header['CRVAL3'])
                    if int(stokes) == 1:
                        stokesid = "i"
                    elif int(stokes) == 2:
                        stokesid = "q"
                    elif int(stokes) == 3:
                        stokesid = "u"
                    elif int(stokes) == 4:
                        stokesid = "v"  
                        
                print(stokesid)  
            
                # Store cutouts to outdir. These are saved automatically.
                fits.writeto('%s/vlass_cutout%s_%s_%s_%s_%s.fits'%(cutout_outdir, epoch, tile[i], 
                    subtile[i], source_pos, stokesid), trimmed.data, trimmed.header, overwrite=True)  
        else:
            print('>>>>> Image/cube tile ID %s and subtile ID %s is not '\
            'available to create a cutout.'%(tile[i], subtile[i]))
            #sys.exit()
            
    if do_mosaic:
        if source_dec < 0:
            source_pos = '%s%s'%(source_ra, source_dec)
        if source_dec > 0:
            source_pos = '%s+%s'%(source_ra, source_dec)
            
        outprefix = 'vlass_cutout_mosaic_%s_%s_%s'%(epoch, source_pos, stokesid)
        mosaic_images(cutouts, outdir=cutout_outdir, outprefix=outprefix)
        
        
 
if __name__ == "__main__":

    descStr = """ VLASS Cutout Tool/code. """
    parser = argparse.ArgumentParser(description=descStr)
    
    parser.add_argument("-e", dest="epoch", type=int, default=2,
                        help="Choose epoch. e.g. 1, 2 or 3.")
    parser.add_argument("-ra", dest="source_ra", type=float,
                        help="Right ascension of the source in degrees.")
    parser.add_argument("-dec", dest="source_dec", type=float,
                        help="Declination of the source in degrees.")   
    parser.add_argument("-s", dest="cut_size", type=float, default=3,
                        help="The size of the cutout in arcminutes. Default is 3.")
    parser.add_argument("-fits_idir", dest="fits_indir", type=str, 
                        help="Give path to input directory where the fits file to"
                        " generate cutouts are stored.")
    parser.add_argument("-fits_odir", dest="fits_outdir", type=str, 
                        help="Give path to a directory to store the final cutouts/mosaiced fits images")
    parser.add_argument("-do_msc", dest="do_mosaic", action='store_true',
                        help="If enabled, a mosaic will be generated. ")     
                                                                         
    args = parser.parse_args()      
    
    
    source_ra = args.source_ra #12.007905 #12.085246752041655 #  12.34 #
    source_dec = args.source_dec #-2.356234 #-2.2197933226664 # -3.015 #
    epoch = args.epoch
    cutout_size = args.cut_size/60. #Cutout size is 3arcmin

    # Read subtile catalog for a specified epoch.  
    cirada_subtile = glob.glob('subtile_catalogs/CIRADA_VLASS%dQ*_subtile.csv'%epoch)[0]
    
    csv_data = read_csv_subtiles(cirada_subtile)

    subtile_crval1 = csv_data['CRVAL1']
    subtile_crval2 = csv_data['CRVAL2']

    #convert to normal array from csv-like array. 
    subtile_crval1 = numpy.asarray([i for i in subtile_crval1])
    subtile_crval2 = numpy.asarray([i for i in subtile_crval2])
 
    # determine the adjacent subtiles. 
    nearby_subtile_idx = find_nearest_subtiles(source_ra, source_dec, subtile_crval1, subtile_crval2)
    nearby_subtile_data = csv_data.iloc[nearby_subtile_idx]
    
    # determine the subtiles for which the source appears. Note that for VLASS, a source 
    # can appear in more than one tile due to overlapping of the edges.
    source_overlap_subtiles = overlap_subtile(source_ra, source_dec, nearby_subtile_data) 
     
    #TODO: we also want to return near by subtiles for which the required cutout sizes overlaps.
    
    cutout_image(source_overlap_subtiles, cutout_size, fits_indir=args.fits_indir, 
       cutout_outdir=args.fits_outdir, do_mosaic=args.do_mosaic)
       










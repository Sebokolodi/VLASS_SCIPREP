#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
According to the user manual at 
https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide,
the central position of each cube (in the CRVAL1 and CRVAL2 keywords in the headers)
differs slightly per channel due to the image-plane correction for w-terms. 
 
The majority of this script is obtained from process_cube.py written by Mark Lucy.
    
Script written by Lerato Baidoo, 09 May 2023.
"""


import glob
import numpy as np
import os
import sys
import argparse

       
def convolve(image, bmaj, bmin, bpa, outfile):

    """ Convolve the input image to the desired resolution. 
    
    image: input image to convolve
    bmaj, bmin, bpa: desired beam major, minor axis and
                     position angle.
    Outfile: name of the convolved casa image file.
    
    """

    # read the beam sizes from the images, to ensure that
    # the values provided are feasible. If bmaj/bmin are smaller
    # than those of the image, we will get an error. 
    hdrdict = imhead(im) 
    try:
        # if it hasn't been convolved before it will have these:   
        # This is only reading for Stokes I. There are values for Stokes Q and U.
        # So far it seems like these are the same. Double check with the VLASS team.
        image_bmaj = hdrdict['perplanebeams']['beams']['*0']['*0']['major']['value']
        image_bmin = hdrdict['perplanebeams']['beams']['*0']['*0']['minor']['value']

    except:
        # if it has been convolved already, but we want to convolve it to our desired beam.
        image_bmaj = hdrdict['restoringbeam']['major']['value']
        image_bmin = hdrdict['restoringbeam']['minor']['value']


    if image_bmaj > bmaj or image_bmin > bmin:
        sys.exit('Specified beam maj, min is smaller (%.2f,%.2f) deg than the image beam (%.2f,%.2f) deg.'\
                %(bmaj, bmin, image_bmaj, image_bmin))

    else:
        print('Convolving the data...')
        imsmooth(imagename=image, major=bmaj, minor=bmin, pa=bpa, 
             targetres=True, outfile=outfile, overwrite=True) 
             
             

def align_images(image, outregrid):

    """ Aligns the images to the corrected point centres (crvals).
    
    image: image to regrid
    outregrid: output name for a regridded image.
    
    """

    hist = imhistory(im) # read history
    sub = 'Uncorrected CRVAL' # get uncorrected crvlist
    crvlist = [s for s in hist if sub in s]
    
    # there are situations where the image way already be realigned. 
    # so the crvlist will be empty.
    if bool(crvlist):
        print('Running a regriding to the corrected pointing centre.')
        crv1 = crvlist[0]
        crv2 = crvlist[1]
    
        orig_crv1 = float(crv1.split(' ')[3])*np.pi/180.
        orig_crv2 = float(crv2.split(' ')[3])*np.pi/180.    

        template = imregrid(im, template='get')
        print(template['csys']['direction0']['crval'])
        template['csys']['direction0']['crval'] = [orig_crv1, orig_crv2]

        # regrid to the new direction.
        imregrid(image, output=outregrid, template=template)
        # export CASA image to fits file.
        exportfits(outregrid, outregrid +'.fits', overwrite=True) 

    else:
        print('No need for regridding')
        exportfits(image, outregrid+'.fits')



if __name__ == "__main__":

    descStr = """
    Pre-processes VLASS data. First, the script convolves to a user desired resolution,  
    then repositions the images to the correct central values. 
    """
    parser = argparse.ArgumentParser(description=descStr)
    parser.add_argument("-indir", dest="input_directory", type=str, default='./',
                        help="Give name of input directory. Include path where necessary.")
    parser.add_argument("-img", dest="imagename", type=str, default=None,
                        help="Image name to process. The image must be inside the indir.")
    parser.add_argument("-outdir", dest="output_directory", type=str, default='./',
                        help="Give name of output directory to save outputs maps. " 
                        "Include path where necessary.")
    parser.add_argument("-bmaj", dest="beam_major_axis", type=float, default=5,
                        help="The desired beam major axis in arcseconds.")
    parser.add_argument("-bmin", dest="beam_minor_axis", type=float, default=5,
                        help="The desired beam minor axis in arcseconds.")
    parser.add_argument("-bpa", dest="beam_posang_axis", type=float, default=0,
                        help="The desired beam position angle in degrees")   
    parser.add_argument("-all", dest="process_all_img", action='store_true',
                        help="If enabled, the code will process all images inside indir.")
                           
                     
    args = parser.parse_args()
    
    majbm = args.beam_major_axis
    minbm = args.beam_minor_axis
    posangle = args.beam_posang_axis
          
    if args.process_all_img:
        imlist = glob.glob( args.input_directory + '/'+ '*.fits')
        for im in imlist:
            outsmooth = im.split('.fits')[0] + '.conv'
            outregrid = outsmooth + '.regrid'
            # the code requires that the images are convolved first before regrided. Otherwise,
            # CASA complains. But what happens when the images are already convolved?
            convolve(image=im, bmaj=majbm, bmin=minbm, bpa=posangle, outfile=outsmooth)
            align_images(outsmooth, outregrid)   
        
             
    else:
        imname =  args.imagename  
        if imname is None:
            sys.exit('No image is provided nor an option to run all images.')
        else:   
            im = args.input_directory + '/' +  imname
        
        outsmooth = im.split('.fits')[0] + '.conv'
        outregrid = outsmooth + '.regrid'
        convolve(image=im, bmaj=majbm, bmin=minbm, bpa=posangle, outfile=outsmooth)
        align_images(outsmooth, outregrid)

    
    
    
    

import glob
import time
from astropy.io import fits 


#Instructions to installing astropy inside CASA. 
#import pip._internal
#pip._internal.main(['install', 'astropy', '--user'])
#import subprocess, sys
#subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'astropy'])
#https://stackoverflow.com/questions/52289107/installation-of-astropy-in-casa (installing pip worked)
#https://astropy-cjhang.readthedocs.io/en/latest/install.html (install astropy worked)


"""
This script converts Stokes cubes into frequency cubes.

The input cubes from the VLASS are Stokes cubes, that is, per 
spectral window, the cube has Stokes I, Q and U.

1. first we concatenate the cubes into Stokes and frequency cube.
2. then we split the mega cube along Stokes parameters: I, Q and U.

This script uses CASA packages.
"""

imlist = glob.glob('J005000-023000/VLASS2.1.cc.T10t02.J005000-023000.06.2048.v1.spw*.IQU.iter3.image.pbcor.tt0.subim.conv.regrid.fits')

def order_images_freq(imlist):
    """Order input cubes by increasing frequencies.
    
    imlist: the list of image cubes to convert into frequency cubes.
    """
    
    ref_freq = []
    for im in imlist:
        hdr = imhead(im)
        ref_freq.append(hdr['refval'][2])

    ordered_imlist = [x for _,x in sorted(zip(ref_freq, imlist))]
    ordered_freq = sorted(ref_freq)
    
    return ordered_freq, ordered_imlist
    
    
def add_freq_hdr(primary_image, frequencies):
    """Add frequencies to the fits header. This is important 
    for non-regular frequencies (e.g. in case of missing spectral windows).
    
    primary_image: mega cube image.
    frequencies: frequencies of the present spectral windows
    """

    #useful link: 
    #https://stackoverflow.com/questions/36458109/save-2-fits-files-with-3-headers-astropy-or-pyfits
    image = fits.open(primary_image)
    primary_data = image[0].data
    primary_hdr = image[0].header
    
    new_hdu = fits.HDUList()
    new_hdu.append(fits.PrimaryHDU(data=primary_data,header=primary_hdr))
    new_hdu.append(fits.ImageHDU(frequencies, name='Frequencies [Hz]'))
    new_hdu.writeto(primary_image, overwrite=True)
    
    return primary_image

    

def concatenate_images(imlist, outfile):    

    """Concatename the Stokes cubes into Stokes and frequency cube.

    imlist: a list of images to convert into frequency cubes.
    outfile: a file to use for concatenated cube (the mega cube)
    """
    #TODO: add frequencies in the header extension file. 
    
    start = time.time()
    freqs, ordered_cubes = order_images_freq(imlist)
    #turn a list of ordered cubes into a string suitable for 
    #https://casa.nrao.edu/docs/casaref/image.imageconcat.html .
    
    infiles = ' '.join(ordered_cubes)
    
    print('>>> Concatenating the Stokes cubes into Stokes-Freq cubes.')
    # then concatenate them
    ia.imageconcat(infiles=infiles, outfile=outfile, relax=True)
    #exportfits(outfile, outfile+'.fits')
    end = time.time()
    print('>>> Contenating completed. Execution time = %.2f secs'%(end-start))  
       
    
    return freqs
 

def convert_to_freqCubes(imlist, outfile):

    
    freqs = concatenate_images(imlist, outfile=outfile) 
    
    print('>>> Splitting mega cube image into Stokes I, Q and U.')
    imsubimage(imagename=outfile, outfile=outfile + '_i', stokes='I')
    imsubimage(imagename=outfile, outfile=outfile + '_q', stokes='Q')
    imsubimage(imagename=outfile, outfile=outfile + '_u', stokes='U')
    
    #check if there are complete number of channels, if some are missing, 
    #replace them with nans.
    
    
    exportfits(outfile + '_i', outfile + '_i.fits')
    exportfits(outfile + '_q', outfile + '_q.fits')
    exportfits(outfile + '_u', outfile + '_u.fits')
    
    #print('>> Adding frequency list to the header extension 1.')
    #add_freq_hdr(outfile + '_i.fits', freqs)
    #add_freq_hdr(outfile + '_q.fits', freqs)
    #add_freq_hdr(outfile + '_u.fits', freqs)
    #print('>>> Done adding frequency list to header fits header file.') 
    
    #TODO: I need to add the frequency list at the final stages, because casa does not propagate the header extension. 
    # Alternatively is to add the frequency info as a comment/history. 
    

    
convert_to_freqCubes(imlist, outfile='image_test')   


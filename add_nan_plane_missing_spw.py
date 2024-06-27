import glob
import numpy 
from astropy.io import fits

#Purpose: to generate nan plane images. 

#frequency_dict = {'spw2':2.028, 'spw3': 2.156, 'spw4': 2.284, 'spw5':2.412, 'spw6':2.54, 'spw7':2.668, 
#        'spw8':2.796, 'spw9':2.924, 'spw10':3.052, 'spw11': 3.18, 'spw12':3.308, 'spw13':3.436, 
#        'spw14':3.564, 'spw15':3.692, 'spw16':3.82, 'spw17':3.948}


def open_image(image):

    """Read the central frequencies"""

    img = fits.open(image)[0]
    header = img.header
    data =  img.data

    return header


def update_header(header, freq0):

    """Update image header"""

    hdr = header.copy()
    hdr['CRVAL3'] = freq0*1e9
    return hdr


def main(image_lists):

    frequency_dict = {'spw2':2.028, 'spw3': 2.156, 'spw4': 2.284, 'spw5':2.412, 'spw6':2.54, 'spw7':2.668, 
        'spw8':2.796, 'spw9':2.924, 'spw10':3.052, 'spw11': 3.18, 'spw12':3.308, 'spw13':3.436, 
        'spw14':3.564, 'spw15':3.692, 'spw16':3.82, 'spw17':3.948}

    central_frequencies = [open_image(x)['CRVAL3']/1e9 for x in image_list]

    missing_freq_dict = {}
    for i in range(2, 18):
        spectral_window = 'spw%d'%i
    
        key = spectral_window
        value = frequency_dict[key]
        if value not in central_frequencies:
            missing_freq_dict.update({key:value})

    
    # if there are missing spectral windows
    if len(missing_freq_dict) > 0:
        print('The missing spectral windows: ')
        print(missing_freq_dict)
        #get the central frequency of the  missing spectral window
        ref_image = image_list[0]
        ref_hdr = open_image(ref_image)
        ref_naxis1 = ref_hdr['naxis1']
        ref_naxis2 = ref_hdr['naxis2']
        ref_naxis3 = ref_hdr['naxis3']
        ref_naxis4 = ref_hdr['naxis4']
        ref_freq0 = round(ref_hdr['crval3']/1e9, 3)

        ref_key = [i for i in frequency_dict if frequency_dict[i]==ref_freq0]
 
        for k in missing_freq_dict.keys():

            outname = image_list[0].replace(ref_key[0], k) 
            freq0 = missing_freq_dict[k]
            new_ref_hdr = update_header(ref_hdr, freq0)
        
            image_cube = numpy.zeros([ref_naxis4, ref_naxis3, ref_naxis1, ref_naxis2], dtype=numpy.float32)    
            image_cube[:, 0, :, :] = float(numpy.nan)
            fits.writeto(outname, image_cube, new_ref_hdr, overwrite=True)
            print('Done creating a missing spectral window')
            
             


path_to_tile = '/home/lerato/vlass/cubes/T10t35/'
#subtiles = ['J005000-023000', 'J005800-013000','J010200-013000' ,'J010201-033000', 'J011001-023000',  'J011800-003000']
#subtiles = ['J184200-033000']
#subtiles = ['J204600-023000', 'J210201-033000']
subtiles = ['J224959-003000', 'J225400-013000','J230600-003000']

for subtile in subtiles:

    image_list = glob.glob(path_to_tile + subtile + '/' + 'conv/' + '*tt0.subim.conv.regrid.fits')
    print(image_list)
    main(image_list)

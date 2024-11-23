from astropy.io import fits

files = ['RT_amax_2cm.fits']


for infits in files:

    # Fits file that works...
    hdu_works = fits.open('header_works.fits')
    hdr_works = hdu_works[0].header
    #print(hdr_works)
    data = hdu_works[0].data
    print(data.shape)

    # Fits file that's wrong... 
    hdu_wrong = fits.open(infits)
    hdr_wrong = hdu_wrong[0].header
    out_data = hdu_wrong[0].data[0,0,:,:]
    print(out_data.shape)

    # Make a new header with structure of works, with numbers from wrong... 
    new_hdr = hdr_works
    new_hdr['NAXIS'] = 2
    new_hdr['NAXIS1'] = hdr_wrong['NAXIS1']
    new_hdr['NAXIS2'] = hdr_wrong['NAXIS2']

    new_hdr['CTYPE1'] = hdr_wrong['CTYPE1']
    new_hdr['CRVAL1'] = hdr_wrong['CRVAL1']
    new_hdr['CRPIX1'] = hdr_wrong['CRPIX1']
    new_hdr['CDELT1'] = hdr_wrong['CDELT1']

    new_hdr['CTYPE2'] = hdr_wrong['CTYPE2']
    new_hdr['CRVAL2'] = hdr_wrong['CRVAL2']
    new_hdr['CRPIX2'] = hdr_wrong['CRPIX2']
    new_hdr['CDELT2'] = hdr_wrong['CDELT2']

    new_hdr['RESTFREQ'] = hdr_wrong['RESTFREQ']

    fits.writeto(infits[:-5]+'_fixed.fits', out_data, new_hdr, overwrite=True)



























'''

fits_file = 'RT.fits'

fits.info(fits_file)



with fits.open(fits_file, 'update') as f:
    for hdu in f:
        hdu.header.remove('FLUX_1')
        hdu.header.remove('FLUX_2')
        hdu.header.remove('FLUX_3')
        hdu.header.remove('FLUX_4')




        hdu.header['OBJECT'] = 'CAT'








hdu.header.remove('CRPIX3')
hdu.header.remove('CRVAL3')
hdu.header.remove('CDELT3')
hdu.header.remove('CUNIT3')
hdu.header.remove('CTYPE3')
'''

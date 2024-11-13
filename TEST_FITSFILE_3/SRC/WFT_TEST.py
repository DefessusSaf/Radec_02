import astropy.io.fits as pyfits
import astropy.units as _u
from os import system
from astropy import coordinates
import TicToc
from os import listdir, getcwd
import numpy as np

DIR_IN = 'LOAD_FILE'
DIR_OUT = 'TMP'
fname = listdir(DIR_IN)[0]

fname_input = f'{DIR_IN}/{fname}'
fname_output = f'{DIR_OUT}/TMP_{fname}'

header = pyfits.getheader(fname_input, ext=0)
image_data = pyfits.getdata(fname_input, ext=0)

PIXSIZE = header['XPIXSZ']
WIDTH = header['NAXIS1']
HEIGHT = header['NAXIS2']
FOCALLEN = 551 #header['FOCALLEN']
BINNING = header['XBINNING']
SCALE = BINNING * PIXSIZE/FOCALLEN*_u.rad.to('arcsec')*_u.um.to('mm')

RA = header['OBJCTRA']
DEC = header['OBJCTDEC']

COORDS = coordinates.SkyCoord(RA, DEC, unit=('hour', 'deg'), frame='icrs', equinox='J2000')

print(f'RA = {RA}, DEC = {DEC}')
print(f'RA(deg.) = {COORDS.ra.deg:.5f},  DEC(deg.) = {COORDS.dec.deg:.5f}')

ASTROMETRY_CONFIG = 'CONFIGS/astrometry.cfg'
DOWNSAMPLE = 4
SEARCH_RAD = SCALE*_u.arcsec.to('deg')*np.max([WIDTH, HEIGHT])*1.2

str1 = f'solve-field --config {ASTROMETRY_CONFIG} --overwrite --downsample {DOWNSAMPLE} --cpulimit 3600 --no-plots '
str2 = f'--scale-units app --scale-low {SCALE*0.95:.2f} --scale-high {SCALE*1.05:.2f} '
str3 = f'--ra {COORDS.ra.deg:.5f} --dec {COORDS.dec.deg:.5f} --radius {SEARCH_RAD} '
str4 = f'--x-column X_CENTER --y-column Y_CENTER --sort-column AREA --width {WIDTH} --height {HEIGHT} '

cmd = str1 + str2 + str3 + str4 + 'TMP/XY.fits'

print(cmd)

try:
    TicToc.tic()
    system(cmd)
    TicToc.toc()

    cmd = f'/usr/local/astrometry/bin/new-wcs -i {fname_input} -w TMP/XY.wcs -o TMP/WCS{fname}'
    system(cmd)

    header = pyfits.getheader(f'TMP/WCS{fname}', ext=0)

    fits = pyfits.open(fname_input)
    fits[0].header = header
    fits[0].data = image_data 
    fits.writeto(fname_output, output_verify='silentfix', overwrite=True)
    fits.close()
except:
    print(f'No WCS solution for {fname_input}')
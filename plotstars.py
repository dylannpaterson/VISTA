from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

import gc

from os import listdir,chdir

import scipy.ndimage as sp

import matplotlib.pyplot as plt

import numpy as np
import numpy.ma.mrecords as mrecords

filelist = listdir('/home/dpaterson/Documents/VISTA/Ks')

chdir('/home/dpaterson/Documents/VISTA/Ks')

ii = 0
for files in filelist:
	
	hdusource = fits.open(files)
	if ii == 0:
		sourcetab = hdusource[1].data
	else:
		sourcetab = np.hstack((sourcetab,hdusource[1].data))
	sourcehdr = hdusource[1].header
	ii += 1

hst = np.histogram2d(sourcetab['GLON'] - 360*(sourcetab['GLON'] > 180),sourcetab['GLAT'],[1200,190])

plt.imshow(np.transpose(hst[0]),origin='lower')
#plt.contour((hst[1][:180]+hst[1][1:])/2.0 ,(hst[2][:180] +hst[2][1:])/2.0  ,np.transpose(hst[0]))

plt.show()

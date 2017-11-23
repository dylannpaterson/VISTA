from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

import scipy.optimize as opt

import os

import gc

import scipy.ndimage as sp

import matplotlib.pyplot as plt

import numpy as np
import numpy.ma.mrecords as mrecords

from joblib import Parallel, delayed

import time

def gauss(x,mu,A,sigma):
	return A*np.exp(-np.divide(np.power(x - mu,2),2.0*sigma**2))

def getTileMags(filename):

	hdusource = fits.open(filename)

	sourcetab = hdusource[1].data
	sourcehdr = hdusource[1].header

	if('RADECSYS' in sourcehdr):
    
    		sourcehdr['RADESYSa'] = sourcehdr['RADECSYS']
    		del sourcehdr['RADECSYS']

	apcor = sourcehdr['APCOR3']
	magzp = sourcehdr['MAGZPT']
	
	if('ESO DET DIT' in sourcehdr):
		exptime = sourcehdr['ESO DET DIT']
	else:
		exptime = hdusource[0].header['ESO DET DIT']
		
	if('ESO DET DIT' in sourcehdr):
		tile = sourcehdr['OBJECT']
	else:
		tile = hdusource[0].header['OBJECT']


	if('ESO TEL AIRM START' in sourcehdr):
		airm = 0.5*(sourcehdr['ESO TEL AIRM START'] + sourcehdr['ESO TEL AIRM END'])
	else:
		airm = 0.5*(hdusource[0].header['ESO TEL AIRM START'] + hdusource[0].header['ESO TEL AIRM END'])


	w = wcs.WCS(naxis=2)
	w.wcs.ctype = [sourcehdr['TCTYP3'],sourcehdr['TCTYP5']]
	w.wcs.crpix = [sourcehdr['TCRPX3'],sourcehdr['TCRPX5']]
	w.wcs.crval = [sourcehdr['TCRVL3'],sourcehdr['TCRVL5']]
	w.wcs.pc = [[sourcehdr['TC3_3'],sourcehdr['TC3_5']],[sourcehdr['TC5_3'],sourcehdr['TC5_5']]]

	sourceclass = np.where(sourcetab['Classification']==-1)[0]

	radec = w.wcs_pix2world(np.column_stack((sourcetab['X_coordinate'][sourceclass],sourcetab['Y_coordinate'][sourceclass])),1)

	radecsys = SkyCoord(ra=radec[:,0]*u.degree, dec=radec[:,1]*u.degree, frame='icrs')

	lb = radecsys.galactic

	mag = magzp - 2.5*np.log10(sourcetab['Aper_flux_3'][sourceclass]/exptime) - 0.05*(airm-1.0)  - apcor

	record = np.rec.array([np.asarray(lb.l.degree),np.asarray(lb.b.degree),mag],names=['GLON','GLAT','MAG'])
	

	hdusource.close()
	print tile
	return record
	
def minFind(index,source,nsource):
	
	if np.min(np.sqrt(np.power(source['GLON'] - nsource['GLON'],2) + np.power(source['GLAT'] - nsource['GLAT'],2))) <= (1.0/(3600.0)):

		ind = np.argmin(np.sqrt(np.power(source['GLON'] - nsource['GLON'],2) + np.power(source['GLAT'] - nsource['GLAT'],2)))
		

		return index, ind
		


def remain():
	
	vvvk = getTileMags('ADP.2014-11-12T16:18:07.297.fits')

	vvvj = getTileMags('ADP.2014-11-12T16:18:45.047.fits')
	
	hdu2mass = fits.open('2mass_gc.fits')
	
	mass2 = hdu2mass[1].data
	
	indices = Parallel(n_jobs=2, verbose=5,backend='threading')(delayed(minFind)(index,vvvj,vsource) for index,vsource in enumerate(vvvk))
	
	indices = np.asarray(list(set(indices)))
	
	indices = indices[indices != np.array(None)]
	
	kind,jind = zip(*indices)
	
	ksources = vvvk[list(kind)]
	
	jsources = vvvj[list(jind)]
	
	record = np.rec.array([(ksources['GLON'] + jsources['GLON'])/2.0,(ksources['GLAT'] + jsources['GLAT'])/2.0 ,jsources['MAG'],ksources['MAG'],ksources['GLON'] - jsources['GLON'],ksources['GLAT'] - jsources['GLAT']],names=['GLON','GLAT','JMAG','KMAG','DEL_GLON','DEL_GLAT'])
	
	bt = fits.BinTableHDU(record)

	bt.writeto('vvv_b320.fit',overwrite=True)
	
def main():
	
	hdu2mass = fits.open('2mass_gc.fits')
	
	mass2 = hdu2mass[1].data
	
	hduvvv = fits.open('vvv_b320.fit')
	
	vvv = hduvvv[1].data

	mass2 = np.copy(mass2[mass2['GLAT'] < np.max(vvv['GLAT'])])
	
	mass2 = np.copy(mass2[mass2['GLAT'] > np.min(vvv['GLAT'])])
	
	mass2 = np.copy(mass2[mass2['GLON'] < np.max(vvv['GLON'])])
	
	mass2 = np.copy(mass2[mass2['GLON'] > np.min(vvv['GLON'])])
	
	#plt.plot(vvv['JMAG'] - vvv['KMAG'],vvv['KMAG'],'.')
	
	indices = Parallel(n_jobs=2, verbose=5, backend = 'threading')(delayed(minFind)(index,mass2,vsource) for index,vsource in enumerate(vvv))
	
	indices = np.asarray(list(set(indices)))
	
	indices = indices[indices != np.array(None)]
	
	vind,mind = zip(*indices)
	
	vsources = vvv[list(vind)]
	
	msources = mass2[list(mind)]
	
	record = np.rec.array([(vsources['GLON'] + msources['GLON'])/2.0,(vsources['GLAT'] + msources['GLAT'])/2.0 ,vsources['JMAG'],msources['JMAG'],vsources['KMAG'],msources['KMAG']],names=['GLON','GLAT','JMAG_VVV','JMAG_2MASS','KMAG_VVV','KMAG_2MASS'])
	
	bt = fits.BinTableHDU(record)

	bt.writeto('vvv_2mass_b320.fit',overwrite=True)
	
def calibrate2mass():
	
	plt.close('all')
	
	hdudata = fits.open('vvv_2mass_b320.fit')
	
	hduvvv = fits.open('vvv_b320.fit')
	
	vvv = hduvvv[1].data
	
	hdu2mass = fits.open('2mass_gc.fits')
	
	mass2 = hdu2mass[1].data
	
	mass2 = mass2[mass2['GLAT'] < np.max(vvv['GLAT'])]
	
	mass2 = mass2[mass2['GLAT'] > np.min(vvv['GLAT'])]
	
	mass2 = mass2[mass2['GLON'] < np.max(vvv['GLON'])]
	
	mass2 = mass2[mass2['GLON'] > np.min(vvv['GLON'])]

	print np.max(vvv['GLON']), np.min(vvv['GLON']), np.max(vvv['GLAT']), np.min(vvv['GLAT'])
	
	data = hdudata[1].data
	
	data = np.copy(data[data['KMAG_VVV'] > 12.0])
	
	data = np.copy(data[data['KMAG_VVV'] < 13.0])	
	
	Kcor = np.mean(data['KMAG_2MASS'] - data['KMAG_VVV'])
	
	Jcor = np.mean(data['JMAG_2MASS'] - data['JMAG_VVV'])
	
	print Kcor, Jcor
	
	vvv['KMAG'] = vvv['KMAG'] + Kcor
	
	vvv['JMAG'] = vvv['JMAG'] + Jcor
	
	vvv = np.rec.array([vvv['GLON'],vvv['GLAT'],vvv['JMAG'],vvv['KMAG']], dtype = [('GLON',float),('GLAT',float),('JMAG',float),('KMAG',float)])
	
	mass2 = np.rec.array([mass2['GLON'],mass2['GLAT'],mass2['JMAG'],mass2['KMAG']],names = ['GLON','GLAT','JMAG','KMAG'], dtype = [('GLON',float),('GLAT',float),('JMAG',float),('KMAG',float)])
	
	sources = np.hstack((vvv[vvv['KMAG'] >= 12.0],mass2[mass2['KMAG'] < 12.0]))

	print len(vvv[vvv['KMAG'] >= 12.0])

	window = sources[np.abs(sources['GLON'] - np.mean(sources['GLON'] - 0.5)) <= 5.0/60.0]

	window = window[np.abs(window['GLAT'] - np.mean(window['GLAT'] - 0.5)) <= 5.0/60.0]
	
	hist,xedges,yedges = np.histogram2d(window['JMAG'] - window['KMAG'], window['KMAG'] ,[100,100],normed=True)
	
	#plt.imshow(np.transpose(hist),origin = 'lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],interpolation='nearest',cmap = 'nipy_spectral')

	#plt.contour((xedges[1:] + xedges[:-1])/2.0,(yedges[1:] + yedges[:-1])/2.0,np.transpose(hist),10,colors = 'k',linewidths=2.0)

	#plt.plot(window['JMAG'] - window['KMAG'], window['KMAG'],'g.')

	#plt.xlim([0.8,1.6])
	#plt.ylim([16.0,11.0])
	#plt.gca().set_aspect('auto')

	#plt.show()

	bt = fits.BinTableHDU(sources)

	bt.writeto('b320.fit',overwrite=True)
	
	#baade = sources[np.abs(sources['GLON'] - 1.14) <= 5.0/60.0]

	#baade = baade[np.abs(baade['GLAT'] + 4.18) <= 5.0/60.0]

	#plt.figure()
	#plt.plot(baade['JMAG'] - baade['KMAG'],baade['KMAG'],'k.')

	#plt.xlim([0.6, 1.4])
	#plt.ylim([16,11])

	#baade = baade[np.logical_and(baade['KMAG'] < 14.0,baade['KMAG'] > 12.4)]

	#baade = baade[np.logical_and(baade['JMAG'] - baade['KMAG'] > 0.775,baade['JMAG'] - baade['KMAG'] < 1.15)]

	#plt.figure()
	#plt.hist(baade['JMAG'] - baade['KMAG'], bins = 15,histtype = 'step',color= 'k')

	#plt.xlim([0.6, 1.4])

	#print np.median(baade['JMAG'] - baade['KMAG'])

	#hst,bins = np.histogram(baade['JMAG'] - baade['KMAG'], 15)

	#bin_centre = (bins[:-1] + bins[1:])/2.0

	#popt,pov = opt.curve_fit(gauss,bin_centre,hst,p0 = [1.0,300.0,0.3])

	#print popt

	#plt.plot(np.linspace(0.6,1.4,100),gauss(np.linspace(0.6,1.4,100),*popt))

	#plt.show()
	
remain()
main()	
calibrate2mass()

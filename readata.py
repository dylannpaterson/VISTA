# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:24:24 2017

@author: dpaterson
"""

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

import os

import gc

import scipy.ndimage as sp

import matplotlib.pyplot as plt

import numpy as np
import numpy.ma.mrecords as mrecords

from joblib import Parallel, delayed

import time


def getTileMags(filename):

	hdusource = fits.open(filename)

#hdusourceJ = fits.open('/home/dpaterson/Downloads/b250/ADP.2017-01-18T11:58:38.802.fits')
	sourcetab = hdusource[1].data
	sourcehdr = hdusource[1].header

	if('RADECSYS' in sourcehdr):
    
    		sourcehdr['RADESYSa'] = sourcehdr['RADECSYS']
    		del sourcehdr['RADECSYS']

	apcor = sourcehdr['APCOR3']
	magzp = sourcehdr['MAGZPT']
	exptime = sourcehdr['ESO DET DIT']
	tile = sourcehdr['OBJECT']

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

	record = np.rec.array([np.asarray(lb.l.degree),np.asarray(lb.b.degree),mag],names=['GLON','GLAT','KMAG'])
	

	hdusource.close()
	print tile
	return(tile,record)

def RmDuplicates(tilelist,recarray):
	
	for ii in range(len(tilelist)):

		print tilelist[ii]

		tileid = int(tilelist[ii][1:])-200

		neighbours = FindNeighbours(tileid)

		for nb in neighbours:

			print '-', nb  
			newnb = FindDuplicates(recarray[ii],recarray[np.where(tilelist==nb)[0][0]],ii+1,int(nb[1:])-200);
			recarray[ii] = newnb
	
		bt = fits.BinTableHDU(recarray[ii])

		bt.writeto(os.getcwd()+'/Ks/'+tilelist[ii]+'.fits',overwrite=True)

	return tilelist,recarray

def minFind(index,source,tileb):

	if np.min(np.sqrt(np.power(source['GLON'] - tileb['GLON'],2) + np.power(source['GLAT'] - tileb['GLAT'],2))) < 1e-4:

		return(index)
			
def FindDuplicates(tilea,tileb,tileida,tileidb):

	dupesind = []
	ii = 0
	jj = 0

	lonlowr,lonupr = np.percentile(tilea['GLON'],[10,90])

	latlowr,latupr = np.percentile(tilea['GLAT'],[15,85])


	if tileidb - tileida == 1:

		mask = tilea['GLON'] > lonupr

	elif tileidb - tileida == 14:

		mask = tilea['GLAT'] > latupr

	elif tileidb - tileida == 15:

		mask = np.logical_and(tilea['GLON'] > lonupr,tilea['GLAT'] > latupr)

	elif tileidb - tileida == 13:

		mask = np.logical_and(tilea['GLON'] > lonupr,tilea['GLAT'] < latlowr)

	st = time.clock()

	ind = Parallel(n_jobs=4)(delayed(minFind)(index,source,tileb) for index,source in enumerate(tilea[mask]))

	print 'time = ', time.clock() - st

	indf = np.asarray(list(set(ind)))

	indf = indf[np.not_equal(indf,None)]

	print 'removed', len(indf), 'duplicates of', np.size(tilea[mask])
	tilearem = np.delete(tilea[mask],indf)

	return(np.hstack((tilea[~mask],tilearem)))
	
def FindNeighbours(tileid):

		neighbours = []

		if tileid >14 and tileid <183:

			neighbours.append('b'+str(tileid+214))
			neighbours.append('b'+str(tileid-14+200))

			if tileid % 14 != 0 and tileid % 14 != 1:

				neighbours.append('b'+str(tileid+201+14))
				neighbours.append('b'+str(tileid-1+200+14))

				neighbours.append('b'+str(tileid+201-14))
				neighbours.append('b'+str(tileid-1+200-14))

				neighbours.append('b'+str(tileid+201))
				neighbours.append('b'+str(tileid-1+200))

			elif tileid % 14 == 0:

				neighbours.append('b'+str(tileid-1+200+14))
				neighbours.append('b'+str(tileid-1+200-14))
				neighbours.append('b'+str(tileid-1+200))

			elif tileid % 14 == 1:

				neighbours.append('b'+str(tileid+201+14))
				neighbours.append('b'+str(tileid+201-14))
				neighbours.append('b'+str(tileid+201))

		elif tileid <= 14:

			neighbours.append('b'+str(tileid+214))

			if tileid % 14 != 0 and tileid % 14 != 1:

				neighbours.append('b'+str(tileid+201+14))
				neighbours.append('b'+str(tileid-1+200+14))

				neighbours.append('b'+str(tileid+201))
				neighbours.append('b'+str(tileid-1+200))

			elif tileid % 14 == 0:

				neighbours.append('b'+str(tileid-1+200+14))
				neighbours.append('b'+str(tileid-1+200))


			elif tileid % 14 == 1:

				neighbours.append('b'+str(tileid+201+14))
				neighbours.append('b'+str(tileid+201))


		elif tileid >=183:

			neighbours.append('b'+str(tileid-14+200))

			if tileid % 14 != 0 and tileid % 14 != 1:

				neighbours.append('b'+str(tileid+201-14))
				neighbours.append('b'+str(tileid-1+200-14))

				neighbours.append('b'+str(tileid+201))
				neighbours.append('b'+str(tileid-1+200))

			elif tileid % 14 == 0:

				neighbours.append('b'+str(tileid-1+200-14))
				neighbours.append('b'+str(tileid-1+200))

			elif tileid % 14 == 1:

				neighbours.append('b'+str(tileid+201-14))
				neighbours.append('b'+str(tileid+201))

		completed = []

		for i in range(1,tileid):

			completed.append('b'+str(200+i))


		return list(set(neighbours) - set(completed))

def main():
	filepath = '/media/dpaterson/Seagate Expansion Drive/Science/vvv/Ks2/'
 
	tilenames = open('Ksfilenames_version3.txt','r')

	tilelist = []
	reclist = []

	for rows in tilenames:
		tilelist.append(rows[:-1])

	ii = 0

	[reclist.append(getTileMags(filepath+tile+'.fits')) for tile in tilelist]

	tileids = np.asarray(zip(*reclist)[0])
	recarrays = list(zip(*reclist)[1])

	tiles,reclistdr = RmDuplicates(tileids,recarrays)
	
	#reclist = list(zip(*reclist)[1])


	ii = 0

	for arr in reclistdr:

		if ii == 0:
			catlg = arr
		else:
			catlg = np.concatenate((catlg,arr))

		ii +=1

	bt = fits.BinTableHDU(catlg)

	bt.writeto('vvv_Ks.fits',overwrite=True)


main()

#sourceJtab = hdusourceJ[1].data
#sourceJhdr = hdusourceJ[1].header


#if('RADECSYS' in sourceJhdr):
    
#    sourceJhdr['RADESYSa'] = sourceJhdr['RADECSYS']
#    del sourceJhdr['RADECSYS']

#apcorJ = sourceJhdr['APCOR3']
#magzpJ = sourceJhdr['MAGZPT']
#exptimeJ = sourceJhdr['ESO DET DIT']
#skynoiseJ = sourceJhdr['SKYNOISE']
#rcoreJ = sourceJhdr['RCORE']

#if('HIERARCH ESO TEL AIRM START' in sourceJhdr):
#	airmJ = 0.5*(sourceJhdr['HIERARCH ESO TEL AIRM START'] + sourceJhdr['HIERARCH ESO TEL AIRM END'])
#else:
#	airmJ = 0.5*(hdusourceJ[0].header['HIERARCH ESO TEL AIRM START'] + hdusourceJ[0].header['HIERARCH ESO TEL AIRM END'])



#wJ = wcs.WCS(naxis=2)

#wJ.wcs.ctype = [sourceJhdr['TCTYP3'],sourceJhdr['TCTYP5']]
#wJ.wcs.crpix = [sourceJhdr['TCRPX3'],sourceJhdr['TCRPX5']]
#wJ.wcs.crval = [sourceJhdr['TCRVL3'],sourceJhdr['TCRVL5']]
#wJ.wcs.pc = [[sourceJhdr['TC3_3'],sourceJhdr['TC3_5']],[sourceJhdr['TC5_3'],sourceJhdr['TC5_5']]]

#classJ = np.where(sourceJtab['Classification']==-1)[0]


#x = np.zeros((len(classKs),2))
#y = np.zeros((len(classKs),2))


#xJ = np.zeros((len(classJ),2))
#yJ = np.zeros((len(classJ),2))



#maglimKs = magzpKs - 2.5*np.log10(5*np.sqrt((np.pi*rcoreKs**2))*skynoiseKs/exptimeKs)- 0.05*(airmKs-1.0)  - apcorKs

#xKs[:,0] = sourceKtab['X_coordinate'][classKs]
#xKs[:,1] = sourceKtab['Y_coordinate'][classKs]
#yKs[:,0] = sourceKtab['X_coordinate'][classKs] + sourceKtab['X_coordinate_err'][classKs]
#yKs[:,1] = sourceKtab['Y_coordinate'][classKs] + sourceKtab['Y_coordinate_err'][classKs]

    
	
#radecerrKs = wKs.wcs_pix2world(yKs,0) - radecKs

#magJ = magzpJ - 2.5*np.log10(sourceJtab['Aper_flux_3'][classJ]/exptimeJ) - 0.05*(airmJ-1.0) - apcorJ

#maglimJ = magzpJ - 2.5*np.log10(5*np.sqrt((np.pi*rcoreJ**2))*skynoiseJ/exptimeJ)- 0.05*(airmJ-1.0)  - apcorJ

#xJ[:,0] = sourceJtab['X_coordinate'][classJ]
#xJ[:,1] = sourceJtab['Y_coordinate'][classJ]
#yJ[:,0] = sourceJtab['X_coordinate'][classJ] + sourceJtab['X_coordinate_err'][classJ]
#yJ[:,1] = sourceJtab['Y_coordinate'][classJ] + sourceJtab['X_coordinate_err'][classJ]
    

    
#radecJ = wJ.wcs_pix2world(xJ,0)
#radecerrJ = wJ.wcs_pix2world(yJ,0) - radecJ

#xypixKs = wKs.wcs_world2pix(radecKs[:],0)-1

#xypixJ = wKs.wcs_world2pix(radecJ[:],0)-1
    
#srclistJ = []
#srclistKs = []
    
#for ii in range(len(classJ)):
    
    #psind = np.where((np.sqrt((np.power(radecKs - radecJ[ii],2))[:,0] + (np.power(radecKs - radecJ[ii],2))[:,1]) < psf))[0]
    
    
    #if (psind.size != 0):

        #srclistJ.append(ii)
        #srclistKs.append(psind[0])
            


#xypixKsr = wKs.wcs_world2pix(radecKs[srclistKs],0)-1

#xypixJr = wKs.wcs_world2pix(radecJ[srclistJ],0)-1



#magJ = magzpJ - 2.5*np.log10(sourceJtab['Aper_flux_3'][classJ][srclistJ]/exptimeJ) - 0.05*(airmJ-1.0) - apcorJ



#hst = np.histogram2d((magJ - magKs),magKs,180,normed=True)

#print np.min(magKs), np.max(magKs), maglimKs
#print np.min(magJ), np.max(magJ),maglimJ
#plt.gca().invert_yaxis()

#plt.show()

#plt.figure()

#plt.contour((hst[1][:180]+hst[1][1:])/2.0 ,(hst[2][:180] +hst[2][1:])/2.0  ,np.transpose(hst[0]),[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

#plt.xlim([0.0,2.0])
#plt.ylim([11,15.5])

#plt.xlabel('J - Ks')
#plt.ylabel('Ks')

#plt.gca().invert_yaxis()

#plt.show()

#plt.figure()

#plt.imshow(hst[0])

#plt.show()
#    for ii in range(len(sourceJtab)):
#        print abs(radecKs[jj][0] - radecJ[ii][0]),radecerrKs[jj][0]+radecerrJ[ii][0], abs(radecJ[ii][1] - radecKs[jj][1]), radecerrKs[jj][1]+radecerrJ[ii][1]
#        
#        if abs(radecKs[jj][0] - radecJ[ii][0]) <= (radecerrKs[jj][0]+radecerrJ[ii][0]): #and abs(radecJ[ii][1] - radecKs[jj][1]) <= (radecerrKs[jj][1]+radecerrJ[ii][1])):
#            srclist.append([magKs[jj],magJ[ii]])
#            print magKs[jj],magJ[ii]
#            radecJ = np.delete(radecJ,(ii),axis=0)
#            break
#    if(len(radecJ) == 0):
#        break
#    
#plt.plot((magJ - magKs),magKs)
#plt.imshow(abs(imagedat[0:100][0:100]),cmap='binary',norm=LogNorm(vmin=5000,vmax=10000))
#plt.plot(xx[:],yy[:],'.')

#del imagedat

#hduimage.close()

#hdusourceJ.close()

#gc.collect()

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:03:39 2017

@author: dpaterson
"""

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

import csv

cfile = open('fp_2mass.fp_psc31933.csv')

creader = csv.reader(cfile,delimiter=' ', skipinitialspace=True)

mass = np.zeros(shape=(3e7,4))

for index,row in enumerate(creader):
	mass[index] = [float(row[0]), float(row[1]), float(row[2]), float(row[4])]

	
mass = mass[:index]
radecsys = SkyCoord(ra=mass[:,0]*u.degree, dec=mass[:,1]*u.degree, frame='icrs')

lb = radecsys.galactic

record = np.rec.array([np.asarray(lb.l.degree),np.asarray(lb.b.degree),mass[:,2],mass[:,3]],names=['GLON','GLAT','JMAG','KMAG'])

record = record[np.logical_or(record['GLON'] > 349.5,record['GLON'] < 10.5)]

record = record[record['GLAT'] > -10.5]

record = record[record['GLAT'] < 10.5]
bt = fits.BinTableHDU(record)

bt.writeto('2mass_gc.fits',overwrite=True)

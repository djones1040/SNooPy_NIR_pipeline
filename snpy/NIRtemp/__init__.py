#!/usr/bin/env python
'''A python module to generate NIR lightcurve templates
from Avelino et al. (2019)
Author:	 David Jones
Version:  0.0

'''
import sys,os,string
import numpy as num
from scipy.interpolate import interp1d
import scipy.optimize
import pickle

debug=0

template_bands = ['Y','J','H','K']

base = os.path.dirname(globals()['__file__'])
if base == '':	base = '.'

def dm152s(dm15):
   '''Convert from dm15 parameter to stretch.'''
   return((3.06-dm15)/2.04)

class dm15_template:
	def __init__(self):
		
		self.dm15 = None
		self.normalize = 1   # Do we force max of lightcurve = 0?
		self.t = None
		self.interp = {}

		self.setup_interpolators()

		self.Y = None;  self.eY = None
		self.J = None;  self.eJ = None
		self.H = None;  self.eH = None
		self.K = None;  self.eK = None
	  
	def __getstate__(self):

		d = self.__dict__.copy()
		d['interp'] = {}
		return d

	def setup_interpolators(self):
		
		x,y,z = num.loadtxt(os.path.join(base, "Y_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['Y'] = interp1d(x,y, bounds_error=False)
		self.interp['eY'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "J_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['J'] = interp1d(x,y, bounds_error=False)
		self.interp['eJ'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "H_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['H'] = interp1d(x,y, bounds_error=False)
		self.interp['eH'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "K_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['K'] = interp1d(x,y, bounds_error=False)
		self.interp['eK'] = interp1d(x,z, bounds_error=False)
		
	def mktemplate(self, dm15, dm15_int=None, dm15_colors='int', generate=0):
		'''2nd and 3rd arguments ignored.'''
		self.dm15 = dm15

		if generate:
			# generate the model light-curve from -10 to 80 in 1 day increments.
			self.t = num.arange(-15,81, 1.0)
			for band in ['Y','J','H','K']:
				self.__dict__[band],self.__dict__['e'+band], mask = \
					self.eval(band, self.t)
			#self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)
			
	def deltaTmax(self, band):
		'''Given the current dm15, what is the time of maximum of [band]
		relative to B-band.'''
		if band == 'Y':  return 0.0
		if band == 'J':  return 0.0
		if band == 'H':  return 0.0
		if band == 'K':  return 0.0
		return 0

	def domain(self, band):
		'''returns the valid domain of the template'''
		if band not in self.interp:
			self.set_interpolators()
			s = dm152s(self.dm15)
		return (self.interp[band].x.min()*s, self.interp[band].x.max()*s)

	def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
		'''Evaluate the template in band [band] at epochs [times].  Optionally
		redshift by (1+[z]).	If [mag]=1, return in magnitudes, otherwise return
		in flux units.  If [sextrap]=1, extrapolate beyond the training sample
		by using a stretch.  Use [gen] to specifiy the generation of the template.
		If you want the Tmax - Tmax(B) offset applied, set [toff] to True,
		otherwise, Tmax will be at 0 for every filter.'''

		if band not in self.interp:
			self.setup_interpolators()

		if toff:
			evt = (times - self.deltaTmax(band))/(1+z)
		else:
			evt = times/(1+z)

		# Apply a stretch consistent with fits to dm15
		s = dm152s(self.dm15)
		evm = self.interp[band](evt/s)
		eevm = self.interp['e%s'%band](evt/s)#num.zeros(evm.shape)   # No uncertainties
		mask = ~num.isnan(evm)
		evm[~mask] = -1.0
		if mag:
			return evm,eevm,mask
		else:
			return num.power(10, -0.4*evm),eevm,mask

class st_template:
	def __init__(self):
		self.st = None
		self.normalize = 1   # Do we force max of lightcurve = 0?
		self.t = None
		self.interp = {}

		self.setup_interpolators()
		
		self.Y = None;  self.eY = None
		self.J = None;  self.eJ = None
		self.H = None;  self.eH = None
		self.K = None;  self.eK = None

	def __getstate__(self):
		
		d = self.__dict__.copy()
		d['interp'] = {}
		return d

	def setup_interpolators(self):
		
		x,y,z = num.loadtxt(os.path.join(base, "Y_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['Y'] = interp1d(x,y, bounds_error=False)
		self.interp['eY'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "J_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['J'] = interp1d(x,y, bounds_error=False)
		self.interp['eJ'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "H_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['H'] = interp1d(x,y, bounds_error=False)
		self.interp['eH'] = interp1d(x,z, bounds_error=False)
		x,y,z = num.loadtxt(os.path.join(base, "K_Template_phase_mu_blue.dat"),
							unpack=True)
		self.interp['K'] = interp1d(x,y, bounds_error=False)
		self.interp['eK'] = interp1d(x,z, bounds_error=False)
		
	def mktemplate(self, st, dm15_int=None, dm15_colors='int', generate=0):
		'''2nd and 3rd arguments ignored.'''
		self.st = st
		
		if generate:
			# generate the model light-curve from -10 to 80 in 1 day increments.
			self.t = num.arange(-15,81, 1.0)
			for band in ['Y','J','H','K']:
				self.__dict__[band],self.__dict__['e'+band], mask = \
					self.eval(band, self.t)
				#self.__dict__['e'+band] = num.where(mask, self.__dict__['e'+band], -1.0)

	def deltaTmax(self, band):
		'''Given the current dm15, what is the time of maximum of [band]
		relative to B-band.'''
		if band == 'Y':  return 0.0
		if band == 'J':  return 0.0
		if band == 'H':  return 0.0
		if band == 'K':  return 0.0
		return 0

	def domain(self, band):
		'''returns the valid domain of the template'''
		s = dm152s(self.dm15)
		if band not in self.interp:
			self.set_interpolators()
		return (self.interp[band].x.min()*s, self.interp[band].x.max()*s)

	def eval(self, band, times, z=0, mag=1, sextrap=1, gen=1, toff=True):
		'''Evaluate the template in band [band] at epochs [times].  Optionally
		redshift by (1+[z]).  If [mag]=1, return in magnitudes, otherwise return
		in flux units.  If [sextrap]=1, extrapolate beyond the training sample
		by using a stretch.  Use [gen] to specifiy the generation of the template.
		If you want the Tmax - Tmax(B) offset applied, set [toff] to True,
		otherwise, Tmax will be at 0 for every filter.'''

		if band not in self.interp:
			self.setup_interpolators()

		if toff:
			evt = (times - self.deltaTmax(band))/(1+z)
		else:
			evt = times/(1+z)

		# Apply a stretch consistent with fits to dm15
		evm = self.interp[band](evt/self.st)
		eevm = self.interp['e%s'%band](evt/self.st)#num.zeros(evm.shape)   # No uncertainties
		mask = ~num.isnan(evm)
		evm[~mask] = -1.0
		if mag:
			return evm,eevm,mask
		else:
			return num.power(10, -0.4*evm),eevm,mask


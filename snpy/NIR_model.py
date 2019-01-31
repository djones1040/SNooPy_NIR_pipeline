'''Model.py:  a module that defines the SN models to be fit by SNOOPY.

A base class (Model) is defined to handle most of the heavy-lifting and boiler
plate around scipy.optimize.leastsq. Defining a new model is done by sub-
classing Model and overriding the member functions.

New:  Add an optional [decline_param] to choose between a dm15 model and stretch
	  (st)	model'''
import os,string
from snpy import ubertemp,NIR_ubertemp
from snpy import kcorr
from snpy.utils import redlaw
from numpy.linalg import cholesky
from scipy import stats
from scipy.optimize import leastsq
from scipy.optimize import brent
import scipy.interpolate
from numpy import *
#from numpy import median, bool, diag
from numpy.linalg import inv
import pickle
from model import model
from model import read_table, base

class EBV_NIR_model(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The parameters you can fit:

   - dm15 (decline rate)
   - Tmax (time of peak B maximum)
   - DM	  (distance modulus)
   - EBVhost  (host galaxy extinction)

   The model is constructed by assuming a peak B absolute magnitude	 and B-X
   colors based on the current value of dm15.  The colors are from Folatelli
   et al. (2010), as are the calibration of Bmax vs dm15.  For the latter,
   there are 6 calibrations, based on the sample used to make the fit.	The
   default is 6 (best observed, excluding heavily extincted SNe), but you can
   choose a different calibration by setting that argument in the fit() call.
   Aside from the instrinsic colors, a global extinction parameter EBVhost
   is applied to each light-curve, as well as Milky way extinction from 
   the SN object's EBVgal.	The value of R_V for the host galaxy is
   not a parameter, but is controled by the choice of calibration in order to
   remain consistent with Folatelli et al. (2009).	The R_V for the galactic
   extinction is taken from the SN object (default 3.1).'''

   def __init__(self, parent, stype='dm15'):

	  if stype != 'dm15':
		 raise ValueError, "This model only supports the dm15 parameter"
	  model.__init__(self, parent)
	  self.rbs = ['u','B','V','g','r','i','Y','J','H','K','Bs','Vs','Rs','Is',
			'J_K','H_K']
	  self.parameters = {'DM':None, 'dm15':None, 'EBVhost':None, 'Tmax':None}
	  self.errors = {'DM':0, 'dm15':0, 'EBVhost':0, 'Tmax':0}
	  self.template = NIR_ubertemp.template()
	  # R_V as a function of which calibration fit number (see Folatelli et
	  #	 al. (2009) table 9
	  self.Rv_host = {1:0, 2:3.10, 3:1.50, 4:1.46, 5:1.46, 6:1.01}
	  self.dRv_host = {1:0, 2:0, 3:0.11, 4:0.10, 5:0.33, 6:0.31}
	  self.M0 = {1:-19.07,2:-19.39,3:-19.15,4:-19.13,5:-19.16,6:-19.11}
	  self.dM0 = {1:0.01, 2:0.02, 3:0.02, 4:0.01, 5:0.03, 6:0.02}
	  self.b = {1:1.03, 2:0.98, 3:0.94, 4:1.05, 5:0.94, 6:1.08}
	  self.db = {1:0.25, 2:0.41, 3:0.11, 4:0.11, 5:0.12, 6:0.11}
	  # B-X pseudo-colors from Folatelli et al. (2009) table 3
	  self.colors = {'u':-0.32, 'B':0.0, 'V':-0.02, 'g':0.05, 'r':-0.09, 'i':-0.63,
			'Y':-0.69, 'J':-0.65, 'H':-0.79, 'K':-0.61, 'J_K':-0.65, 'H_K':-0.79}
	  self.dcolors = {'u':0.04, 'B':0, 'V':0.01, 'g':0.02, 'r':0.02, 'i':0.02,
			'Y':0.03, 'J':0.02, 'H':0.03, 'K':0.05, 'J_K':0.02, 'H_K':0.03,
			'Bs':0, 'Vs':0, 'Rs':0, 'Is':0}
	  self.color_slopes = {'u':-0.47, 'B':0.0, 'V':0.12, 'g':0.05, 'r':0.29,
						   'i':0.39, 'Y':0.63, 'J':0.67, 'H':0.66, 'K':0.26,
						   'J_K':0.67, 'H_K':0.66}
	  self.dcolor_slopes = {'u':0.25, 'B':0, 'V':0.05, 'g':0.06, 'r':0.07,
							'i':0.08, 'Y':0.17, 'J':0.10, 'H':0.11, 'K':0.18,
							'J_K':0.10, 'H_K':0.11, 'Bs':0, 'Vs':0, 'Rs':0,
							'Is':0}
	  self.do_Robs = 0
	  self.Robs = {}


   def setup(self):
	  if 'EBVhost' not in self.args:
		 if len(self._fbands) < 2:
			raise RuntimeError, "Error:	 to solve for EBVhost, you need to fit more than one filter"

	  self.calibration = self.args.get('calibration',6)
	  self.gen = self.args.get('gen',2)
	  for band in self._fbands:
		 #cal = self.args.get('cal',6)
		 cal = self.calibration
		 self.Robs[band] = kcorr.R_obs(band, self.parent.z, 0, 0.01, 0,
			   self.Rv_host[cal], self.parent.Rv_gal, self.parent.k_version,
			   redlaw=self.parent.redlaw)
	  
   def guess(self, param):
	  s = self.parent
	  if param == 'Tmax':
		 Tmaxs = []
		 for f in s.data:
			Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
		 return median(Tmaxs)

	  if param == 'DM':
		 # Quick DM based on Ho = 72
		 if s.z < 1e-10:
			raise ValueError, "SN redshift is too close to zero.  Set it properly."
		 return 43.11 + 5*log10(s.z)

	  if param == 'dm15':
		 # choose just the average dm15:
		 return(1.1)

	  return(0.0)

   def __call__(self, band, t, extrap=False):
	  self.template.mktemplate(self.dm15)
	  if len(shape(t)) == 0:
		 t = array([t])
	  t = t - self.Tmax
	  rband = self.parent.restbands[band]

	  # Now build the lc model
	  temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
			gen=self.gen, extrap=extrap)
	  K,mask2 = self.kcorr(band, t)
	  temp = temp + K

	  # Apply reddening correction:
	  # Figure out the reddening law
	  if self.do_Robs:
		 self.Robs[band] = kcorr.R_obs(band, self.parent.z, t, self.EBVhost,
			   self.parent.EBVgal, self.Rv_host[self.calibration], 
			   self.parent.Rv_gal, self.parent.k_version, 
			   redlaw=self.parent.redlaw)
		 temp = temp + self.Robs[band]*(self.EBVhost + self.parent.EBVgal)
	  else:
		 # Apply Robs*EBVgal:
		 R = self.MWR(band, t)
		 temp = temp + self.Robs[band]*self.EBVhost + R*self.parent.EBVgal
	  temp = temp + self.DM + self.MMax(rband, self.calibration)

	  return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
	  Tmaxs = []
	  Mmaxs = []
	  eMmaxs = []
	  rbands = []
	  self.template.mktemplate(self.dm15)
	  for band in bands:
		 rband = self.parent.restbands[band]
		 # find where the template truly peaks:
		 x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
		 Tmaxs.append(x0 + self.Tmax)
		 mmax = self.DM + self.MMax(rband, self.calibration)
		 if not restframe and band in self.parent.ks_tck:
			# add the K-correction
			mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
		 if not deredden:
			if self.do_Robs:
			   Robs = kcorr.R_obs(band, self.parent.z, x0, self.EBVhost, 
					 self.parent.EBVgal, self.Rv_host[self.calibration], 
					 self.parent.Rv_gal, 'H3', redlaw=self.parent.redlaw)
			   mmax = mmax + Robs*(self.EBVhost + self.parent.EBVgal)
			else:
			   # Apply Robs*EBVgal:
			   if band in self.parent.Robs:
				  if type(self.parent.Robs[band]) is type(()):
					 R = scipy.interpolate.splev(Tmaxs[-1], self.parent.Robs[band])
				  else:
					 R = self.parent.Robs[band]
			   else:
				  EBV = max(self.parent.EBVgal, 0.01)
				  R = kcorr.R_obs(band, self.parent.z, 0, 0, EBV,
								  self.Rv_host[self.calibration], 
								  self.parent.Rv_gal,
								  self.parent.k_version,
								  redlaw=self.parent.redlaw)
			   EBV = max(self.EBVhost, 0.01)
			   Robs = kcorr.R_obs(band, self.parent.z, 0, EBV,
					 0, self.Rv_host[self.calibration], self.parent.Rv_gal,
					 self.parent.k_version, redlaw=self.parent.redlaw)
			   mmax = mmax + Robs*self.EBVhost + R*self.parent.EBVgal
		 Mmaxs.append(mmax)
		 eMmaxs.append(self.errors['DM'])
		 rbands.append(rband)
	  return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def MMax(self, band, calibration=6):
	  '''Given self.dm15, return the absolute magnitude at maximum for the given
	  filter [band].  The calibration paramter allows you to choose which
	  fit (1-6) in Folatelli et al. (2009), table 9'''
	  if band == 'Bs':
		 return -19.319 + (self.dm15-1.1)*0.634
	  elif band == 'Vs':
		 return -19.246 + (self.dm15-1.1)*0.606
	  elif band == 'Rs':
		 return -19.248 + (self.dm15-1.1)*0.566
	  elif band == 'Is':
		 return -18.981 + (self.dm15-1.1)*0.524
	  elif band in ['u','B','V','g','r','i','Y','J','H','K','J_K','H_K']:
		 return self.M0[calibration] + (self.dm15-1.1)*self.b[calibration] -\
				self.colors[band] - self.color_slopes[band]*(self.dm15 -1.1)
	  else:
		 return -19.0

   def systematics(self, calibration=6, include_Ho=False):
	  '''Returns the systematic errors in the paramters as a dictionary.  
	  If no estimate is available, return None for that paramter.'''
	  systs = dict.fromkeys(self.parameters.keys())
	  # DM contains systematics for Ho, plus average of calibration
	  #	 uncertainties
	  syst_DM = []
	  syst_EBV = []
	  weights = []
	  for band in self._fbands:
		 rb = self.parent.restbands[band]
		 # make a call to get the weights
		 mod,err,mask = self.__call__(band, self.parent.data[band].MJD)
		 weights.append(sum(where(mask, power(err,-2), 0)))
		 ddm15 = self.dm15 - 1.1
		 Robs = kcorr.R_obs(band, self.parent.z, 0, self.EBVhost,
			   self.parent.EBVgal, self.Rv_host[calibration], 
			   self.parent.Rv_gal, self.parent.k_version, 
			   redlaw=self.parent.redlaw)
		 dRobs = Robs*self.dRv_host[calibration]/self.Rv_host[calibration]
		 syst_DM.append(power(ddm15*self.db[calibration],2)+\
						power(ddm15*self.dcolor_slopes[rb],2)+\
						power(self.EBVhost*dRobs, 2)+\
						power(self.dM0[calibration], 2)+\
						power(self.dcolors[rb]*Robs, 2) +\
						power(0.06*Robs, 2) +\
						#power(2.17*velerr/(3e5*self.parent.z),2) +\
						power(0.06,2))
		 syst_EBV.append(power(0.06*Robs,2))
	  syst_DM = array(syst_DM)
	  weights = array(weights)
	  systs['DM'] = sum(weights*syst_DM)/sum(weights)
	  if include_Ho:
		 systs['DM'] += power(2.17*0.1,2)	  # assume 10% error in Ho
	  systs['DM'] = sqrt(systs['DM'])
	  systs['EBVhost'] = 0.06
	  systs['dm15'] = 0.06
	  systs['Tmax'] = 0.34
	  return(systs)

class EBV_NIR_model2(model):
   '''This model fits any number of lightcurves with CSP uBVgriYJHK templates
   or Prieto BsVsRsIs templates.  The parameters you can fit:

   - dm15 or st (decline rate or stretch)
   - Tmax (time of peak B maximum)
   - DM	  (distance modulus)
   - EBVhost  (host galaxy extinction)

   The model is constructed by assuming a peak absolute magnitudes in tall
   the filters, as derived in Burns et al. 2011.  Calibrations were determined
   using MCMC modeling on all filters at once, determining M_X and b_X for
   each filter, and one value for R_V.	There are 2 calibrations, based on the 
   sample used to make the fit and the prior used on the extinction.  Default
   is 0, where the 2 red SNe are excluded and the blue sub-sample is used to
   anchor the colors.  Value of 1 is for the sample where the two red
   SNe were included.  A global extinction parameter EBVhost
   is applied to each light-curve, as well as Milky way extinction from 
   the SN object's EBVgal.	The value of R_V for the host galaxy is
   not a parameter, but is controled by the choice of calibration
   R_V for the galactic extinction is taken from the SN object (default 3.1).'''

   def __init__(self, parent, stype='st'):

	  if stype not in ['dm15','st']:
		 raise ValueError, "This model only supports the dm15 and st parameter"
	  model.__init__(self, parent)
	  self.rbs = ['u','B','V','g','r','i','Y','J','H']
	  self.parameters = {'DM':None, stype:None, 'EBVhost':None, 'Tmax':None}
	  self.errors = {'DM':0, stype:0, 'EBVhost':0, 'Tmax':0}
	  if stype == 'dm15':
		 self.template = NIR_ubertemp.template()
	  else:
		 self.template = NIR_ubertemp.stemplate()
	  self.stype = stype
	  
	  if stype in ['st']:
		 self.a,self.ea,self.b,self.eb,self.c,self.ec,self.Rv_host, self.eRv_host,self.sigSN = read_table(os.path.join(base,'st_calibration2.dat'))
	  else:
		 self.a,self.ea,self.b,self.eb,self.c,self.ec,self.Rv_host, self.eRv_host,self.sigSN = read_table(os.path.join(base,'dm15_calibration2.dat'))

	  self.do_Robs = 0
	  self.Robs = {}


   def setup(self):
	  # check to see if we have more than one filter when solving for EBV
	  if 'EBVhost' not in self.args:
		 if len(self._fbands) < 2:
			raise RuntimeError, "Error:	 to solve for EBVhost, you need to fit more than one filter"

	  self.calibration = self.args.get('calibration',0)
	  self.gen = 2

	  for band in self._fbands:
		 self.Robs[band] = kcorr.R_obs(band, self.parent.z, 0, 0.01, 0,
			   self.Rv_host[self.calibration], self.parent.Rv_gal, 
			   self.parent.k_version, redlaw=self.parent.redlaw)
	  
   def guess(self, param):
	  s = self.parent
	  if param == 'Tmax':
		 Tmaxs = []
		 for f in s.data:
			Tmaxs.append(s.data[f].MJD[argmin(s.data[f].mag)])
		 return median(Tmaxs)

	  if param == 'DM':
		 # Quick DM based on Ho = 72
		 if s.z < 0.0001:
			return 0
		 return 43.11 + 5*log10(s.z)

	  if param == 'dm15':
		 # choose just the average dm15:
		 return(1.1)

	  if param == 'st':
		 return (1.0)

	  return(0.0)

   def __call__(self, band, t, extrap=False):
	  self.template.mktemplate(self.parameters[self.stype])
	  t = t - self.Tmax
	  rband = self.parent.restbands[band]

	  # Now build the lc model
	  temp,etemp,mask = self.template.eval(rband, t, self.parent.z, 
			gen=self.gen, extrap=extrap)
	  K,mask2 = self.kcorr(band, t)
	  temp = temp + K

	  # Apply reddening correction:
	  # Figure out the reddening law
	  if self.do_Robs:
		 self.Robs[band] = kcorr.R_obs(band, self.parent.z, t, self.EBVhost,
			   self.parent.EBVgal, self.Rv_host[self.calibration], 
			   self.parent.Rv_gal, self.parent.k_version,
			   redlaw=self.parent.redlaw)
		 temp = temp + self.Robs[band]*(self.EBVhost + self.parent.EBVgal)
	  else:
		 # Apply Robs*EBVgal:
		 R = self.MWR(band, t)
		 temp = temp + self.Robs[band]*self.EBVhost + R*self.parent.EBVgal
	  temp = temp + self.DM + self.MMax(rband, self.calibration)

	  return temp,etemp,mask*mask2

   def get_max(self, bands, restframe=0, deredden=0):
	  Tmaxs = []
	  Mmaxs = []
	  eMmaxs = []
	  rbands = []
	  self.template.mktemplate(self.parameters[self.stype])
	  for band in bands:
		 rband = self.parent.restbands[band]
		 # find where the template truly peaks:
		 x0 = brent(lambda x: self.template.eval(rband, x, gen=self.gen)[0], brack=(0.,5.))
		 Tmaxs.append(x0 + self.Tmax)
		 mmax = self.DM + self.MMax(rband, self.calibration)
		 if not restframe and band in self.parent.ks_tck:
			# add the K-correction
			mmax = mmax + scipy.interpolate.splev(Tmaxs[-1], self.parent.ks_tck[band])
		 if not deredden:
			if self.do_Robs:
			   Robs = kcorr.R_obs(band, self.parent.z, x0, self.EBVhost, 
					 self.parent.EBVgal, self.Rv_host[self.calibration], 
					 self.parent.Rv_gal, 'H3', redlaw=self.parent.redlaw)
			   mmax = mmax + Robs*(self.EBVhost + self.parent.EBVgal)
			else:
			   # Apply Robs*EBVgal:
			   if band in self.parent.Robs:
				  if type(self.parent.Robs[band]) is type(()):
					 R = scipy.interpolate.splev(Tmaxs[-1], self.parent.Robs[band])
				  else:
					 R = self.parent.Robs[band]
			   else:
				  EBV = max(self.parent.EBVgal, 0.01)
				  R = kcorr.R_obs(band, self.parent.z, 0, 0, EBV,
								  self.Rv_host[self.calibration], 
								  self.parent.Rv_gal, self.parent.k_version,
								  redlaw=self.parent.redlaw)
			   EBV = max(self.EBVhost, 0.01)
			   Robs = kcorr.R_obs(band, self.parent.z, 0, EBV,
					 0, self.Rv_host[self.calibration], self.parent.Rv_gal,
					 self.parent.k_version, redlaw=self.parent.redlaw)
			   mmax = mmax + Robs*self.EBVhost + R*self.parent.EBVgal
		 Mmaxs.append(mmax)
		 eMmaxs.append(self.errors['DM'])
		 rbands.append(rband)
	  return(Tmaxs, Mmaxs, eMmaxs, rbands)

   def MMax(self, band, calibration=1):
	  '''Given self.dm15, return the absolute magnitude at maximum for the given
	  filter [band].  The calibration paramter allows you to choose which
	  fit (1-6) in Folatelli et al. (2009), table 9'''
	  if band not in self.a[calibration]:
		 raise ValueError, "Error, filter %s cannot be fit with this calibration"
	  if self.stype in ['st']:
		 delta = self.st - 1.0
	  else:
		 delta = self.dm15 - 1.1
	  return self.a[calibration][band] + self.b[calibration][band]*delta +\
			 self.c[calibration][band]*delta**2

   def systematics(self, calibration=1, include_Ho=False):
	  '''Returns the systematic errors in the paramters as a dictionary.  
	  If no estimate is available, return None for that paramter.'''
	  systs = dict.fromkeys(self.parameters.keys())
	  # DM contains systematics for Ho, plus average of calibration
	  #	 uncertainties
	  syst_DM = []
	  weights = []
	  for band in self._fbands:
		 rb = self.parent.restbands[band]
		 # make a call to get the weights
		 mod,err,mask = self.__call__(band, self.parent.data[band].MJD)
		 weights.append(sum(where(mask, power(err,-2), 0)))
		 if self.stype in ['st']:
			dst = self.st - 1.0
		 else:
			dst = self.dm15 - 1.1
		 Robs = kcorr.R_obs(band, self.parent.z, 0, self.EBVhost,
			   self.parent.EBVgal, self.Rv_host[calibration], 
			   self.parent.Rv_gal, self.parent.k_version,
			   redlaw=self.parent.redlaw)
		 dRobs = Robs*self.eRv_host[calibration]/self.Rv_host[calibration]
		 syst_DM.append(power(dst*self.eb[calibration][rb],2)+\
						power(self.EBVhost*dRobs, 2)+\
						power(self.ea[calibration][rb], 2)+\
						#power(2.17*velerr/(3e5*self.parent.z),2) +\
						power(self.sigSN[calibration][rb],2))
	  syst_DM = array(syst_DM)
	  weights = array(weights)
	  systs['DM'] = sum(weights*syst_DM)/sum(weights)
	  if include_Ho:
		 systs['DM'] += power(2.17*0.1,2)	  # assume 10% error in Ho
	  systs['DM'] = sqrt(systs['DM'])
	  systs['EBVhost'] = 0.06
	  if self.stype == 'dm15':
		 systs['dm15'] = 0.06
	  else:
		 systs['st'] = 0.03
	  systs['Tmax'] = 0.34
	  return(systs)
  

# Fit with Snoopy the light-curve data for a bunch of SNe
# It works for low-z SNe Ia and RAISIN.
#
#--------------------------------------------------------60
# Script created by = 'Arturo_Avelino'
# On date = 2017-01-30 (yyyy-mm-dd)
# last_update = '2019.01.16' # (yyyy-mm-dd)
# version = '0.1.15'
#--------------------------------------------------------60
#
#	USE
#
# i) Open this file in a Text Editor. Select the settings for the fitting.
#
# In the terminal, activate the conda python enviroment where snoopy
# is installed, go to the directory containing this python script and type:
#
#	 python	 ThisScript.py 'dir/where/snoopy/files/are/located'
#
# where the path is the location of the folder containing the snoopy files.
# NOTE: Do not put a slash ("/") at the end of the path.
#
# The output files will be saved in a subdirectory folder called OpticalNIR/Fit/
# that will be created automatically.
# It will be created also a log file listing all the SNe with errors during
# the fitting process.

#--------------------------------------------------------60

import sys # To read arguments in command line
from snpy import *
import numpy as np
import glob # To read the files in my directory
import os # To use command line like instructions
from matplotlib import pyplot as plt
import argparse

#########################################################60
		
class snoopy_fit:
	def __init__(self):
		self.warnings = []

	def addwarning(self,warning):
		print(warning)
		self.warnings.append(warning)

	def check_inputs(self):
		fit_typelist = ('optical','nir','opticalnir')
		if self.options.fit_type not in fit_typelist:
			raise RuntimeError("--fit_type must be in %s"%fit_typelist)
		
	def add_options(self, parser=None, usage=None):
		if parser == None:
			parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

		# The basics
		parser.add_argument('-v', '--verbose', action="count", dest="verbose",
							default=1,help='verbosity level')
		parser.add_argument('--debug', default=False, action="store_true",
							help='debug mode: more output and debug files')
		parser.add_argument('--clobber', default=False, action="store_true",
							help='clobber')

		#BandType
		parser.add_argument('--fit_type', default="opticalnir", type=str,
							help='optical, nir, or opticalnir (default=%default)')
		parser.add_argument('-o','--optical_alone', default=False, action="store_true",
							help=' fit NIR alone, w/o optical (default=%default)')
		parser.add_argument('-k','--no_kcor', default=False, action="store_true",
							help='fit NIR alone, w/o optical (default=%default)')
		parser.add_argument('-k','--no_kcor_stretch', default=False, action="store_true",
							help='don\'t stretch the SED before k-corrections (default=%default)')
		parser.add_argument('--filepath', default=None, type=str,
							help='SNooPy file path, e.g. snoopy_lc/*dat (default=%default)')
		parser.add_argument('--outdir', default='output', type=str,
							help='output directory for fits (default=%default)')
		parser.add_argument('--bandlist', default="gDEC,rDEC,zDEC,iDEC,ps1_g,ps1_r,ps1_i,ps1_z,f125w,f160w", type=str,
							help='list of comma-separated bands to fit, can be empty to fit all bands (default=%default)')

		return parser

	def main(self):

		#- Reading the LC data file names with the snoopy format.
		the_list = glob.glob(self.options.filepath)

		print("# %s SNe in the initial list to be fitted."%len(the_list))

		DirSaveOutput = self.options.outdir #'%s/%s/Fit/'%(
			#os.path.dirname(self.options.filepath),self.options.fit_type)

		#- "If the subdirectory does not exist then create it"
		if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

		#------- Loop over all data -------

		countSN = 0 # Counter number of SNe fitted correctly
		countSNFail = 0 # # Counter number of SNe failed during the fitting.

		#------------------------------

		for file in the_list:
			try:
				print(" ")
				print("\n==================== %s ===================\n"%file[0:14])
				print("%s"%file)
				
				s = get_sn(file)
				s.summary()
				s.choose_model('EBV_NIR_model2')
				
				#- Creation of an array with the name of the filters of this SN.
				FilterNames_array = []
				for band in s.restbands:
					FilterNames_array += [band]
					
				#- Creation of an array with the specific band names to fit for this SN:
				if self.options.bandlist:
					BandsToFit = []; BandsExcludedOfFit = []
					for band in list(s.data.keys()):
						if band in self.options.bandlist.split(','): BandsToFit += [band]
						else: BandsExcludedOfFit += [band]

				else:
					#- Creation of an array with the name the OPTICAL only and
					# NIR only filters.
					# Observer-frame NIR bands in Andy's compilation or RAISINs
					# (but using Snoopy names).
					All_NIR_bands = ['Y','Ydw','JANDI','J2m', 'J', 'Jrc1', 'Jrc2','Jdw',
									 'HANDI', 'H2m', 'H', 'Hdw', 'KANDI', 'Ks2m', 'K',
									 'f125w', 'f160w']
					OpticalBands = [] # List to put the optical-only bands
					NIRbands = [] # List to put the NIR-only bands
					for band in list(s.data.keys()):
						if band not in All_NIR_bands: OpticalBands += [band]
						else: NIRbands += [band]

				#--- Find out if filters (Bs, Vs) or (B,V) are present in the photometry.
				# If so, do a quick fit (without k-corr) to find T_Bmax,
				# if (Bs, Vs) or (B,V) are -not- present, then fit all the LC
				# data with "s.fit()" and write the name of the SN file in the
				# failure text file.
				#-------
				if ('Bs' and 'Vs') in FilterNames_array:
					# print ("Bands ('Bs','Vs') in: %s \n" % (s.name ))
					# To quickly find the TB_max. It is needed 2 bands necessarily
					s.fit(['Bs','Vs'], dokcorr=0)
				#-------
				elif ('B' and 'V')	in FilterNames_array:
					# print ("Bands ('B','V') in: %s \n" % (s.name ))
					# To quickly find the TB_max. It is needed 2 bands necessarily
					s.fit(['B','V'], dokcorr=0)
					#-------
				elif ('B' and 'V0')	 in FilterNames_array:
					# print ("Bands ('B','V') in: %s \n" % (s.name ))
					# To quickly find the TB_max. It is needed 2 bands necessarily
					# s.fit(['B','V0'], dokcorr=0)
					s.fit(['B','V0'])
					#-------
				elif ('B' and 'V1')	 in FilterNames_array:
					# print ("Bands ('B','V') in: %s \n" % (s.name ))
					# To quickly find the TB_max. It is needed 2 bands necessarily
					s.fit(['B','V1'], dokcorr=0)
					#-------
				elif ('BANDI' and 'VANDI')	in FilterNames_array:
					# print ("Bands ('BANDI','VANDI') in: %s \n" % (s.name ))
					# To quickly find the TB_max. It is needed 2 bands necessarily
					s.fit(['BANDI','VANDI'], dokcorr=0)
					#-------
				else: # When there is NOT the bands (B, V) nor (Bs, Vs) in the LC data
					# print ("No bands in %s \n" % (s.name ))
					print("No B,V observed-frame bands found, but don't worry.")

				print("%s. Prefitting (B,V) with no kcorrections: done."%s.name)

				#-------------------------------------------------------------------
				# MAIN FITTING, either, Optical only, Optica+NIR, or specific bands only:

				if self.options.bandlist: # Final fit all the data
					if self.debug: print("Bands to plot:", BandsToFit)
					s.fit(BandsToFit, dokcorr=(not self.options.no_kcor),
						  k_stretch=(not self.options.no_kcor_stretch),
						  reset_kcorrs=True)
					if self.debug: print('Fitted bands:', BandsToFit)

				else:

					if self.options.fit_type == "optical": # Final fit all the data
						if self.debug: print("Bands to plot:", OpticalBands)
						s.fit(OpticalBands, dokcorr=(not self.options.no_kcor),
							  k_stretch=(not self.options.no_kcor_stretch),
							  reset_kcorrs=True)
						if self.debug: print("Fitting optical bands only: done")

					elif self.options.fit_type == "opticalnir": # Final fit all the data
						if self.debug: print("Bands to plot:", list(s.data.keys()))
						s.fit(dokcorr=(not self.options.no_kcor),
							  k_stretch=(not self.options.no_kcor_stretch),
							  reset_kcorrs=True)
						if self.debug: print("Fitting optical+nir bands: done")

					else:
						raise RuntimeError("NIR alone not yet implemented")

				if self.debug:
					print("%s. Fitting all the bands: done with no issues."%s.name)
				#-------------------------------------------------------------------
				#  Saving the snpy data

				# Removing the words "_snoopy.dat" at the end of the name for each SN.
				NameDataFileToSave = file.split('/')[-1].split('.')[0]

				s.save('%s/%s_1stFit.snpy'%(DirSaveOutput,NameDataFileToSave))
				if self.debug:
					print("%s. The '_1stFit.snpy' file created and saved."%s.name)

				#-------------------------------------------------------------------

				#		PLOTTING

				if self.debug: print("%s. Preparing to plot the fit."%s.name)

				plt.close() # Close any possible plot unfinished/leftover.

				# Plot filters
				s.plot_filters(fill=True, outfile="%s/%s_Filters.png"%(DirSaveOutput,NameDataFileToSave))
				plt.close()
				if self.debug: print("%s. Plot filters: done."%s.name)

				# Plot kcorrs
				s.plot_kcorrs(outfile='%s/%s_PlotKcorrs.png'%(DirSaveOutput,NameDataFileToSave))
				plt.close()
				if self.debug: print("%s. Plot kcorrs: done."%s.name)

				countSN = countSN + 1

				# --- Plotting	--->

				plt.figure()

				x_loc = 5; # days

				#- Determine the y location for the text info. It is going to be below
				#  from the maximum of the last filter to be plotted
				if self.options.bandlist:
					y_loc = s.get_max(bands=BandsToFit[0])[1]+0.6
					x_loc = s.get_max(bands=BandsToFit[0])[0]
				else:
					if self.fit_type == "optical":
						y_loc = s.get_max(bands=s.filter_order[-1])[1]+0.6
						x_loc = s.get_max(bands=s.filter_order[-1])[0]
					else:
						y_loc = s.get_max(bands=s.filter_order[-(len(NIRbands)+1)])[1]+0.6
						x_loc = s.get_max(bands=s.filter_order[-(len(NIRbands)+1)])[0]

				s.plot(epoch=True,xrange=(x_loc-20,x_loc+50),yrange=(y_loc+3,y_loc-2))

				print(x_loc,y_loc)
				try:
					plt.text(-20,y_loc+3,r"""$\Delta$m15 = %.2f $\pm$ %.2f
$z_{\rm hel}$ = %.3f
$\mu$ = %.3f $\pm$ %.3f
$T_{\rm Bmax}$ = %.2f $\pm$ %.3f
E(B-V)$_{\rm MW}$ = %.3f
E(B-V)$_{\rm host}$ = %.3f $\pm$ %.3f"""%(
	s.dm15,s.e_dm15,s.z,s.DM,s.e_DM,s.Tmax,
	s.e_Tmax,s.EBVgal,s.EBVhost,s.e_EBVhost))
				except:
					plt.text(-20,y_loc+3,r"""$s$ = %.2f $\pm$ %.2f
$z_{\rm hel}$ = %.3f
$\mu$ = %.3f $\pm$ %.3f
$T_{\rm Bmax}$ = %.2f $\pm$ %.3f
E(B-V)$_{\rm MW}$ = %.3f
E(B-V)$_{\rm host}$ = %.3f $\pm$ %.3f"""%(
	s.st,s.e_st,s.z,s.DM,s.e_DM,s.Tmax,
	s.e_Tmax,s.EBVgal,s.EBVhost,s.e_EBVhost))
					
				plt.savefig("%s/%s_PlotFitText.png"%(DirSaveOutput,NameDataFileToSave),
							format='png')
				plt.close()
				
				if self.debug: print("%s. Plot fit with text: done."%s.name)
				# <--- Plotting	 ---

				#-----------------------------------------------------------------------

				print('%s: All done with no issues.'%s.name)
			except Exception as err:
				countSNFail = countSNFail + 1
				#Removing the words "snoopy.dat" at the end of the name for each SN:
				NameDataFileToSave = file.split('.')[0]
				print "%s. Failed in some part during running this code."%file
				print "%s"%err
		return countSN, countSNFail
				
if __name__ == "__main__":
	usagestring = """python snoopy_fit.py <options>
example: python snoopy_fit.py --filepath '../snoopy_lc_optnir/*dat' --outdir output_opticalnir'
"""

	snpy = snoopy_fit()
	parser = snpy.add_options(usage=usagestring)
	options = parser.parse_args()

	snpy.options = options
	snpy.debug = options.debug
	
	snpy.check_inputs()
	countSN, countSNFail = snpy.main()

	if len(snpy.warnings):
		print('There were warnings!!')
		print((snpy.warnings))
	
	print("\n\n# -- %i SNe fitted and %i failed --"%(countSN, countSNFail))
		

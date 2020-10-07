import logging
import numpy as np
from astropy.io import fits
import owncloud
import getpass
#import psycopg2 as mdb
#import csv
#import os
#import fnmatch
#import cPickle as pickle
#from scipy import spatial
#import ftputil
#import getpass
#import copy
#from sclip.sclip import sclip
#from scipy.interpolate import UnivariateSpline
import scipy
from scipy import signal
#from pyflann import *
from matplotlib import *
from pylab import *

class spectrum:
	def __init__(self, name, kind='norm', extension=1, wavelength='default', linearize=True):
		
		#set atributes
		if isinstance(name, basestring):
			pass
		else:
			name=str(int(name))
		self.name = name
		self.ccd=int(self.name[-1])
		self.date=int(self.name[:6])
		self.run=int(self.name[6:10])
		self.combine_method=int(self.name[10:12])
		self.pivot=int(self.name[12:15])

		#read spectrum
		self.l=-1
		self.f=-1
		self.fe=-1

		if setup.folder_is_root:
			if name[10:12]=='00':
				#reading uncombined spectra
				path=setup.folder+str(self.date)+'/spectra/all/'+self.name+'.fits'
			elif name[10:12]=='01':
				#readingh combined spectra
				path=setup.folder+str(self.date)+'/spectra/com/'+self.name+'.fits'
			elif name[10:12]=='02':
				#reading spectra combined over several epochs
				path=setup.folder+str(self.date)+'/spectra/com2/'+self.name+'.fits'
			else:
				path=None
		else:
			#reding spectra from an exact folder
			path=setup.folder+self.name+'.fits'
		
		try:#read if it exists
			hdulist = fits.open(path)
		except:#otherwise download
			logging.error('Spectrum %s not found. Enable download to get it from the ftp site.' % self.name)
			raise

		#set l, f, and fe
		# dict converting extension names to extension numbers
		instance={'norm':1, 'normalized':1, 'normalised':1, 'flux':0, 'fluxed':0}

		self.f=hdulist[instance[kind]].data
		self.fe=hdulist[2].data
		self.map_f=hdulist[7].data
		crval=hdulist[instance[kind]].header['CRVAL1']
		crdel=hdulist[instance[kind]].header['CDELT1']
		self.l=np.linspace(crval, crval+crdel*len(self.f), len(self.f))
		self.map_l=np.linspace(crval, crval+crdel*len(self.map_f), len(self.map_f))

		#set radial velocity. This is needed if spectra is required in some other velocity frame than default
		rv=hdulist[0].header['RV']
		if rv=='None':
			self.rv=None
		else:
			self.rv=rv

		#set v_bary. This is needed if spectra is required in some other velocity frame than default
		bary=hdulist[0].header['BARYEFF']
		if bary=='None':
			self.vb=None
		else:
			self.vb=bary
		
		#shift into correct velocity frame
		if wavelength=='default':
			#default velocity frame is as is
			pass
		if wavelength=='bary':
			#barycentric velocity frame is default
			pass
		elif wavelength=='observer':
			#observer velocity frame must be corrected back for v_bary
			if self.vb==None:
				logging.warning('Spectrum %s has no barycentric velocity in the header. Assumibg v_bary=0 km/s' % self.name)
			else:
				self.l=self.l*(1+self.vb/299792.458)
		elif wavelength=='object':
			#object velocity frame must be corrected for radial velocity
			if self.rv==None:
				logging.warning('Spectrum %s has no radial velocity in the header. Assumibg v_r=0 km/s' % self.name)
			else:
				self.l=self.l*(1-self.rv/299792.458)
		else:
			pass

		#linearize
		if linearize==True:
			self.linearize()
		else:
			pass
	
	def linearize(self):
		"""
		take whatever the sampling is and linearize it
		"""
		self.f=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.f)
		self.fe=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.fe)
		self.l=np.linspace(self.l[0],self.l[-1],num=len(self.l))
	
	def logarize(self):
		"""
		take whatever the sampling is and make a logaritmic sampling
		"""
		self.f=np.interp(np.logspace(np.log10(self.l[0]),np.log10(self.l[-1]),num=int(len(self.l)*1.15)),self.l,self.f)
		self.fe=np.interp(np.logspace(np.log10(self.l[0]),np.log10(self.l[-1]),num=int(len(self.l)*1.15)),self.l,self.fe)
		self.l=np.logspace(np.log10(self.l[0]),np.log10(self.l[-1]),num=int(len(self.l)*1.15))
	
	def shift(self, rv, linearize=True):
		"""
		shift the spectrum for radial velocity rv, given in km/s
		if linearize=True, the returned wavelength array will be linearized and flux interpolated
		"""
		l=self.l*(1+rv/299792.458)
		if linearize==True:
			self.f=np.interp(np.linspace(l[0],l[-1],num=len(l)),l,self.f)
			self.fe=np.interp(np.linspace(l[0],l[-1],num=len(l)),l,self.fe)
			self.l=np.linspace(l[0],l[-1],num=len(l))
		else:
			self.l=l

		return self

	def add_noise(self, snr, target_snr,skip=True):
		"""
		Adds poissonian noise to make a spectrum with snr into a spectrum with target_snr.
		"""
		if skip and target_snr>snr:#do not raise error if target_snr is larger than snr. Do nothing instead
			return self

		if target_snr>snr:
			raise RuntimeError('Target SNR cannot be larger than the original SNR.')
		elif target_snr==snr: 
			return self
		else:
			sigma=np.sqrt((1.0/target_snr)**2-(1.0/snr)**2)
			noise=np.random.poisson((1.0/sigma)**2, size=len(self.f))
			noise=noise/((1.0/sigma)**2)
			self.f+=noise
			self.f=self.f-1.0
			self.fe+=1.0/noise
			self.fe=self.fe-1.0

			return self

	def interpolate(self,space):
		"""
		interpolate the spectrum to a wavelength space defined in space. 
		"""
		space=np.array(space)
		self.f=np.interp(space,self.l,self.f)
		self.fe=np.interp(space,self.l,self.fe)
		self.l=space
		return self

	def normalize(self,deg,n,func,sl,su,grow=0,smooth=4e6):
		"""
		calculate the normalization for the spectrum and normalize
		"""
		functions.deg=deg

		if func=='cheb' or func=='chebyshev':
			result=sclip((self.l,self.f),chebyshev,int(n),su=su,sl=sl,min_data=100,verbose=False)
			self.f=self.f/result[0]

		if func=='poly':
			result=sclip((self.l,self.f),poly,int(n),su=su,sl=sl,min_data=100,verbose=False)
			self.f=self.f/result[0]

		if func=='spline':
			functions.smooth=smooth
			result=sclip((self.l,self.f),spline,int(n),su=su,sl=sl,min_data=100,verbose=False)
			self.f=self.f/result[0]

		return result[0]

	def convolve(self,fwhm,extend=False):
		"""
		decrease resolution by convolving the spectrum with a gaussian
		returns the kernel
		"""

		#check if wavelength calibration is linear:
		if (self.l[1]-self.l[0])==(self.l[-1]-self.l[-2]):
			linear=True
		else:
			l=self.l
			self.linearize()
			linear=False

		step=self.l[1]-self.l[0]
		kernel=gauss_kern(fwhm/step)
		add_dim=len(kernel)

		if extend==True:
			f=np.insert(self.f,0,np.ones(add_dim)*self.f[0])
			f=np.append(f,np.ones(add_dim)*self.f[-1])
			fe=np.insert(self.fe,0,np.ones(add_dim)*self.fe[0])
			fe=np.append(fe,np.ones(add_dim)*self.fe[-1])
			self.f=signal.fftconvolve(f,kernel,mode='same')[add_dim:-add_dim]
			self.fe=signal.fftconvolve(fe**2,kernel,mode='same')[add_dim:-add_dim]
			self.fe=np.sqrt(self.fe)

		else:
			self.f=signal.fftconvolve(self.f,kernel,mode='same')
			self.fe=signal.fftconvolve(self.fe**2,kernel,mode='same')
			self.fe=np.sqrt(self.fe)

		if linear==False:
			self.interpolate(l)

		return kernel

	def res_degradation(self,r,target_r):
		"""
		degradate resolution from resolving power r to resolving power target_r
		"""

		if r==target_r:
			pass
		elif target_r>r:
			raise RuntimeError('Cannot increase the resolution.')
		else:
			l=np.average(self.l)
			s=l/r
			s_target=l/target_r

			s_conv=np.sqrt(s_target**2-s**2)

			self.convolve(s_conv,extend=True)

	def equalize_resolution(self):
		"""
		convolves a spectrum with a kernel with a variable width. Works by warping the data, performing the convolution and unwarping the data, so it is vectorized (mostly) and fast
		"""

		#check if wavelength calibration is linear:
		if (self.l[1]-self.l[0])==(self.l[-1]-self.l[-2]):
			linear=True
		else:
			l=self.l
			self.linearize()
			linear=False

		#current sampling:
		sampl=self.l[1]-self.l[0]

		#target sigma coresponds to the R=22000. We want the constant sigma, not constant R, so take the sigma that corresponds to the average R=22000
		s_target=np.ones(len(self.map_l))*np.average(self.map_l)/15000.

		#the sigma of the kernel is:
		s=np.sqrt(s_target**2-self.map_f**2)

		#fit it with the polynomial, so we have a function instead of sampled values:
		map_fit=np.poly1d(np.polyfit(self.map_l, s/sampl, deg=6))

		#create an array with new sampling. The first point is the same as in the spectrum:
		l_new=[self.map_l[0]]

		#and the sigma in pixels with which to convolve is
		sampl=self.l[1]-self.l[0]
		s_new=s/sampl/min(s/sampl)

		#now create the non linear sampling. It should be finer where sigma is larger and sparser where sigma is smaller
		#IS THERE A WAY TO VECTORIZE THIS???
		for i in range(int(len(self.l)*np.max(s_new)*min(s/sampl))):
			l_new.append(l_new[i]+sampl/map_fit(l_new[i]))

		#interpolate the spectrum to the new sampling:
		new_f=np.interp(np.array(l_new),self.l,self.f)
		new_fe=np.interp(np.array(l_new),self.l,self.fe)

		plot(l_new, l_new-np.roll(l_new, 1), 'r,')
		show()

		#the kernel is min(s_orig/sampl) pixels large, so we oversample as little as possible (for the very best precision we should oversample more, but this takes time)
		kernel=gauss_kern(min(s/sampl))

		#convolve the warped spectrum:
		con_f=signal.fftconvolve(new_f,kernel,mode='same')
		con_fe=signal.fftconvolve(new_fe**2,kernel,mode='same')

		#inverse the warping:
		self.f=np.interp(self.l,np.array(l_new),con_f)
		self.fe=np.interp(self.l,np.array(l_new),con_fe)
		self.fe=np.sqrt(self.fe)

		if linear==False:
			self.interpolate(l)

		return self

	def median_filter(self,size, extend=False):
		"""
		do a standard median filtering. give size in Angstroms, will be translated to nearest odd number of pixels.
		"""
		#check if wavelength calibration is linear:
		if (self.l[1]-self.l[0])==(self.l[-1]-self.l[-2]):
			linear=True
		else:
			linear=False
			l=self.l
			self.linearize()

		step=self.l[1]-self.l[0]
		add_dim=int(np.ceil(size/step // 2 * 2 + 1))

		if extend==True:
			f=np.insert(self.f,0,np.ones(add_dim)*self.f[0])
			f=np.append(f,np.ones(add_dim)*self.f[-1])
			fe=np.insert(self.fe,0,np.ones(add_dim)*self.fe[0])
			fe=np.append(fe,np.ones(add_dim)*self.fe[-1])
			self.f=signal.medfilt(f,add_dim)[add_dim:-add_dim]
			self.fe=signal.medfilt(fe,add_dim,)[add_dim:-add_dim]
			self.fe=self.fe/np.sqrt(add_dim)

		else:
			self.f=signal.medfilt(self.f,add_dim)
			self.fe=signal.medfilt(self.fe,add_dim)
			self.fe=self.fe/np.sqrt(add_dim)
		
		if linear==False:
			self.interpolate(l)

	def save_fits(self, fname=None):
		"""
		save the spectrum into a 2D fits file
		"""
		if fname==None:
			fname=setup.folder+self.name+'.fits'
		hdu = pyfits.PrimaryHDU(self.f)
		hdu.writeto(fname)
		hdulist = pyfits.open(fname,mode='update')
		prihdr = hdulist[0].header
		prihdr['COMMENT']='File written by galah_tools.py'
		prihdr.set('CRVAL1', self.l[0])
		prihdr.set('CDELT1', self.l[1]-self.l[0])
		prihdr.set('CRPIX1', 1)
		prihdr.set('CUNIT1', 'Angstroms')
		hdulist.flush()
		hdulist.close()
		pyfits.append(fname,self.fe)
		hdulist = pyfits.open(fname,mode='update')
		prihdr = hdulist[1].header
		prihdr['COMMENT']='File written by galah_tools.py'
		prihdr.set('CRVAL1', self.l[0])
		prihdr.set('CDELT1', self.l[1]-self.l[0])
		prihdr.set('CRPIX1', 1)
		prihdr.set('CUNIT1', 'Angstroms')
		hdulist.flush()
		hdulist.close()


	def save_ascii(self, fname=None):
		"""
		save the spectrum into an ascii text file with three columns; wavelength, flux, error
		"""
		if fname==None:
			fname=setup.folder+self.name+'.txt'
		np.savetxt(fname,zip(self.l,self.f,self.fe))

class window_function:
	window=[]
	l=[]
	def __init__(self,l,windows):
		if windows=='':
			window_function.window=np.array([1]*len(l))#if no file is given window function is constant 1
		else:
			window_function.window=read_windows(windows,l)

		window_function.l=l

	def plot_window(self):
		"""
		Plot the window function
		"""
		import matplotlib.pyplot as pl
		fig=pl.figure('Window function (close to continue)')
		pl.plot(window_function.l, window_function.window,'k-')
		pl.xlabel('Wavelength / Angstroms')
		pl.ylabel('Window value')
		pl.title('Window function')
		pl.show()

	def clear(self):
		"""
		Clears a saved window function
		"""
		window_function.window=[]
		window_function.l=[]

class functions:
	def __init__(self):
		deg=0
		smooth=1000

class setup:
	folder=''
	folder_is_root=False
	download=False

	def __init__(self, **kwargs):
		for key in kwargs:
			if key=='folder':
				setup.folder=kwargs['folder']
				if setup.folder[-1]!='/':
					setup.folder=setup.folder+'/'
			elif key=='root_folder':
				setup.folder=kwargs['root_folder']
				setup.folder_is_root=True
				if setup.folder[-1]!='/':
					setup.folder=setup.folder+'/'

			if key=='download':
				setup.download=kwargs['download']


def read_windows(filename,l):
	windows=[]
	cwindows=[]

	#initially the window is set to constant 1
	window=np.ones(len(l))

	#read the file and save windows separately to those specified by lower and upepr limit and those specified by center and width
	writeto=0
	for line in open(filename, 'r'):
		if len(line[:-1])==0: 
			writeto=1
			continue
		if line[0]=='#': pass
		else:
			data=line[:-1].split('\t')
			sign=data[-1][0]
			data=map(float,data)
			data=data+[sign]
			if writeto==0:
				windows.append(data)
			else:
				cwindows.append(data)

	#transform one format into the other
	for i in windows:
		cwindows.append([(i[0]+i[1])/2.0, i[1]-i[0], i[2], i[3], i[4]])

	for w in cwindows:
		if w[4]!='-':
			for n,i in enumerate(l):
				if abs(w[0]-i)<=(w[1]/2.0*(1-w[2])): window[n]*=w[3]
				elif (w[0]+w[1]/2.0*(1-w[2]))<i<(w[0]+w[1]/2.0): window[n]*=(2-2*w[3])/(w[1]*w[2])*i+1-(1-w[3])/w[2]*(2*w[0]/w[1]+1)
				elif (w[0]-w[1]/2.0*(1-w[2]))>i>(w[0]-w[1]/2.0): window[n]*=(2*w[3]-2)/(w[1]*w[2])*i+1-(w[3]-1)/w[2]*(2*w[0]/w[1]-1)
				else: pass
		else:
			for n,i in enumerate(l):
				if abs(w[0]-i)>=(w[1]/2.0): window[n]*=abs(w[3])
				elif (w[0]+w[1]/2.0*(1-w[2]))<i<(w[0]+w[1]/2.0): window[n]*=(2*abs(w[3])-2)/(w[1]*w[2])*i+abs(w[3])-(abs(w[3])-1)/w[2]*(1+2.0*w[0]/w[1])
				elif (w[0]-w[1]/2.0*(1-w[2]))>i>(w[0]-w[1]/2.0): window[n]*=(2-2*abs(w[3]))/(w[1]*w[2])*i+abs(w[3])-(1-abs(w[3]))/w[2]*(2.0*w[0]/w[1]-1)
				else: pass

	return window

def read(name, **kwargs):
	"""
	reads one spectrum
	"""
	spec=spectrum(name, **kwargs)
	return spec

def cc(i1,i2):
	"""
	Cross correlate spectra in array s1 with spectrum s2.
	Returns the shift between spectra (an array if s1 is an array) in km/s
	"""
	s1=copy.deepcopy(i1)
	s2=copy.deepcopy(i2)
	s2.logarize()
	if isinstance(s1, spectrum): s1=[s1]
	ccs=[]
	for s in s1:
		s.interpolate(s2.l)
		ccf=np.correlate(s.f[10:-10]-np.average(s.f[10:-10]),s2.f[10:-10]-np.average(s2.f[10:-10]),mode='same')
		max_index=np.argmax(ccf)
		max_fit=np.polyfit(range(max_index-3,max_index+4),ccf[max_index-3:max_index+4],2)
		max_fit= -max_fit[1]/(max_fit[0]*2.0)
		diff=max_fit-len(ccf)/2.0
		dl=s.l[1]-s.l[0]
		ccs.append(dl*diff/s.l[0]*299792.458)
	
	if len(ccs)==1:
		return ccs[0]
	else:	
		return np.array(ccs)
		
		

def chebyshev(p,ye,mask):
	coef=np.polynomial.chebyshev.chebfit(p[0][mask], p[1][mask], functions.deg)
	cont=np.polynomial.chebyshev.chebval(p[0],coef)
	return cont

def poly(p,ye,mask):
	r=np.polyfit(p[0][mask], p[1][mask], deg=functions.deg)
	f=np.poly1d(r)
	return f(p[0])

def spline(p,ye,mask):
	spl = UnivariateSpline(p[0][mask], p[1][mask],k=functions.deg)
	spl.set_smoothing_factor(5000000)
	return spl(p[0])

def gauss_kern(fwhm):
	""" Returns a normalized 1D gauss kernel array for convolutions """
	size=2*(fwhm/2.355)**2
	size_grid = int(size) # we limit the size of kernel, so it is as small as possible (or minimal size) for faster calculations
	if size_grid<7: size_grid=7
	x= scipy.mgrid[-size_grid:size_grid+1]
	g = scipy.exp(-(x**2/float(size)))
	return g / np.sum(g)

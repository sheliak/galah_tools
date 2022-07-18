import logging
import numpy as np
from astropy.io import fits
#import owncloud
#import getpass
import scipy
from scipy import signal
from matplotlib import *
from pylab import *
import copy

class spectrum:
	def __init__(self, name, kind='norm', extension=1, wavelength='default', linearize=True):
		
		#set atributes
		if isinstance(name, str):
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
			self.hdulist = fits.open(path)
		except:#otherwise download
			logging.error('Spectrum %s not found. Enable download to get it from the ftp site.' % self.name)
			raise

		#set l, f, and fe
		# dict converting extension names to extension numbers
		instance={'norm':1, 'normalized':1, 'normalised':1, 'flux':0, 'fluxed':0}

		self.f=self.hdulist[instance[kind]].data
		self.fe=self.hdulist[2].data
		self.res_map=self.hdulist[7].data
		crval=self.get_param('CRVAL1', ext=instance[kind])
		crdel=self.get_param('CDELT1', ext=instance[kind])
		self.l=np.linspace(crval, crval+crdel*len(self.f), len(self.f))
		self.map_l=np.linspace(crval, crval+crdel*len(self.res_map), len(self.res_map))

		res_b=self.get_param('B', ext=7)
		if res_b is None:
			self.res_b = np.nan
		else:
			self.res_b = float(res_b)

		if len(self.res_map) == 1 or not np.isfinite(self.res_b):
			logging.warning('Resolution changes can not be performed for spectrum %s' % self.name)

		#set radial velocity. This is needed if spectra is required in some other velocity frame than default
		rv=self.get_param('RVCOM')
		if rv=='None':
			self.rv=None
		else:
			self.rv=rv

		#set v_bary. This is needed if spectra is required in some other velocity frame than default
		bary=self.get_param('BARYEFF')
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

	def equalize_resolution(self, r='lowest', precision=0.005, kernel='galah'):
		"""
		Convolves a spectrum with a kernel with a variable width. Works by warping the data, performing the convolution and unwarping the data, so it is vectorized (mostly) and fast

		Parameters:
			r (str or float): resolution of resulting spectra. If r=='lowest', the resulting spectra will have the lowest resolution given in the resolution map
			precision (float): precision of intermitten sampling. Lower number means better precision.
			kernel (str): 'gauss' or 'galah' (default). Will either use a gaussian kernel or a kernel derived for GALAH as the konvolution kernel.
		"""

		#check if wavelength calibration is linear:
		if (self.l[1]-self.l[0])==(self.l[-1]-self.l[-2]):
			linear=True
			l=self.l
		else:
			l=self.l
			self.linearize()
			linear=False

		#current sampling:
		sampl=self.l[1]-self.l[0]

		#target sigma coresponds to the R=22000. We want the constant sigma, not constant R, so take the sigma that corresponds to the average R=22000
		#s_target=np.ones(len(self.map_l))*np.average(self.map_l)/18400.
		if r=='lowest':
			s_target=np.ones(len(self.map_l))*max(self.res_map)*(1+precision)
		elif isinstance(r, (int, long, float)):
			s_target=self.res_map*(1+precision)
			s_target[r>s_target]=r
		else:
			logging.error('Parameter r must be either \'lowest\' or a number.')

		#original sigma
		s_original=self.res_map

		#the sigma of the kernel is:
		s=np.sqrt(s_target**2-s_original**2)


		#fit it with the polynomial, so we have a function instead of sampled values:
		map_fit=np.poly1d(np.polyfit(self.map_l, s, deg=6))

		#create an array with new sampling. The first point is the same as in the spectrum:
		l_new=[self.map_l[0]]

		#minimal needed sampling
		min_sampl=max(s_original)/sampl/sampl

		#keep adding samples until end of the wavelength range is reached
		while l_new[-1]<l[-1]+sampl:
			l_new.append(l_new[-1]+map_fit(l_new[-1])/sampl/min_sampl)

		#interpolate the spectrum to the new sampling:
		new_f=np.interp(np.array(l_new),self.l,self.f)
		new_fe=np.interp(np.array(l_new),self.l,self.fe)

		#plot(l_new, l_new-np.roll(l_new, 1), 'r,')
		#show()

		#convolve
		if kernel=='gauss':
			kernel_=gauss_kern(max(s_original)/sampl)
		elif kernel=='galah':
			kernel_=galah_kern(max(s_original)/sampl, self.res_b)
		else:
			logging.error('Kernel %s is not available. Use gauss or galah.' % kernel)
			return self
		con_f=signal.fftconvolve(new_f,kernel_,mode='same')
		con_fe=signal.fftconvolve(new_fe**2,kernel_,mode='same')

		#inverse the warping:
		self.f=np.interp(self.l,np.array(l_new),con_f)
		self.fe=np.interp(self.l,np.array(l_new),con_fe)
		self.fe=np.sqrt(self.fe)

		if linear==False:
			self.interpolate(l)

		return self

	def synth_resolution_degradation(self, synth, synth_res=300000.0):
		"""
		Take a synthetic spectrum with a very high  resolution and degrade its resolution to the resolution profile of the observed spectrum. The synthetic spectrum should not be undersampled, or the result of the convolution might be wrong.

		Parameters:
			synth np array or similar: an array representing the synthetic spectrum. Must have size m x 2. First column is the wavelength array, second column is the flux array. Resolution of the synthetic spectrum must be constant and higher than that of the observed spectrum.
			synth_res (float): resolving power of the synthetic spectrum

		Returns:
			Convolved syntehtic spectrum as a np array of size m x 2.
		"""

		synth=np.array(synth)
		l_original=synth[:,0]
		#check if the shape of the synthetic spectrum is correct
		if synth.shape[1]!=2: logging.error('Syntehtic spectrum must have shape m x 2.')

		#check if the resolving power is high enough
		sigma_synth=synth[:,0]/synth_res
		if max(sigma_synth)>=min(self.res_map)*0.95: logging.error('Resolving power of the synthetic spectrum must be higher.')

		#check if wavelength calibration of the synthetic spectrum is linear:
		if not (synth[:,0][1]-synth[:,0][0])==(synth[:,0][-1]-synth[:,0][-2]):
			logging.error('Synthetic spectrum must have linear (equidistant) sampling.')		

		#current sampling:
		sampl=synth[:,0][1]-synth[:,0][0]
		galah_sampl=self.l[1]-self.l[0]

		#original sigma
		s_original=sigma_synth

		#required sigma (resample the resolution map into the wavelength range of the synthetic spectrum)
		s_out=np.interp(synth[:,0], self.l, self.res_map)

		#the sigma of the kernel is:
		s=np.sqrt(s_out**2-s_original**2)

		#fit it with the polynomial, so we have a function instead of sampled values:
		map_fit=np.poly1d(np.polyfit(synth[:,0], s, deg=6))

		#create an array with new sampling. The first point is the same as in the spectrum:
		l_new=[synth[:,0][0]]

		#oversampling. If synthetic spectrum sampling is much finer than the size of the kernel, the code would work, but would return badly sampled spectrum. this is because from here on the needed sampling is measured in units of sigma.
		oversample=galah_sampl/sampl*5.0

		#minimal needed sampling
		min_sampl=max(s_original)/sampl/sampl*oversample

		#keep adding samples until end of the wavelength range is reached
		while l_new[-1]<synth[:,0][-1]+sampl:
			l_new.append(l_new[-1]+map_fit(l_new[-1])/sampl/min_sampl)

		#interpolate the spectrum to the new sampling:
		new_f=np.interp(np.array(l_new),synth[:,0],synth[:,1])

		kernel_=galah_kern(max(s_original)/sampl*oversample, self.res_b)
	
		con_f=signal.fftconvolve(new_f,kernel_,mode='same')

		#inverse the warping:
		synth[:,1]=np.interp(l_original,np.array(l_new),con_f)
		return synth


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

	def get_param(self, param, ext=0):
		"""
		Return spectral information or parameter that is written in header. All extensions should have the same set of information.
		"""
		param_get = param.upper()
		if param_get in self.hdulist[ext].header.keys():
			return self.hdulist[ext].header[param_get]
		else:
			logging.error('Parameter %s was not found in extension %d of spectrum %s.' % (param_get, ext, self.name))
			return None

	def close(self):
		"""
		close opened fits file
		"""
		self.hdulist.close()
		return None

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


def sclip(p,fit,n,ye=[],sl=99999,su=99999,min=0,max=0,min_data=1,grow=0,global_mask=None,verbose=True):
		"""
		p: array of coordinate vectors. Last line in the array must be values that are fitted. The rest are coordinates.
		fit: name of the fitting function. It must have arguments x,y,ye,and mask and return an array of values of the fitted function at coordinates x
		n: number of iterations
		ye: array of errors for each point
		sl: lower limit in sigma units
		su: upper limit in sigma units
		min: number or fraction of rejected points below the fitted curve
		max: number or fraction of rejected points above the fitted curve
		min_data: minimal number of points that can still be used to make a constrained fit
		global_mask: if initial mask is given it will be used throughout the whole fitting process, but the final fit will be evaluated also in the masked points
		grow: number of points to reject around the rejected point.
		verbose: print the results or not
		"""

		nv,dim=np.shape(p)

		#if error vector is not given, assume errors are equal to 0:
		if ye==[]: ye=np.zeros(dim)
		#if a single number is given for y errors, assume it means the same error is for all points:
		if isinstance(ye, (int, long, float)): ye=np.ones(dim)*ye
		
		if global_mask==None: global_mask=np.ones(dim, dtype=bool)
		else: pass
		
		f_initial=fit(p,ye,global_mask)
		s_initial=np.std(p[-1]-f_initial)

		f=f_initial
		s=s_initial

		tmp_results=[]

		b_old=np.ones(dim, dtype=bool)

		for step in range(n):
			#check that only sigmas or only min/max are given:
			if (sl!=99999 or su!=99999) and (min!=0 or max!=0):
				raise RuntimeError('Sigmas and min/max are given. Only one can be used.')

			#if sigmas are given:
			if sl!=99999 or su!=99999:
				b=np.zeros(dim, dtype=bool)
				if sl>=99999 and su!=sl: sl=su#check if only one is given. In this case set the other to the same value
				if su>=99999 and sl!=su: su=sl

				good_values=np.where(((f-p[-1])<(sl*(s+ye))) & ((f-p[-1])>-(su*(s+ye))))#find points that pass the sigma test
				b[good_values]=True

			#if min/max are given
			if min!=0 or max!=0:
				b=np.ones(dim, dtype=bool)
				if min<1: min=dim*min#detect if min is in number of points or percentage
				if max<1: max=dim*max#detect if max is in number of points or percentage

				bad_values=np.concatenate(((p[-1]-f).argsort()[-int(max):], (p[-1]-f).argsort()[:int(min)]))
				b[bad_values]=False

			#check the grow parameter:
			if grow>=1 and nv==2:
				b_grown=np.ones(dim, dtype=bool)
				for ind,val in enumerate(b):
					if val==False:
						ind_l=ind-int(grow)
						ind_u=ind+int(grow)+1
						if ind_l<0: ind_l=0
						b_grown[ind_l:ind_u]=False

				b=b_grown

			tmp_results.append(f)

			#check that the minimal number of good points is not too low:
			if len(b[b])<min_data:
				step=step-1
				b=b_old
				break

			#check if the new b is the same as old one and break if yes:
			if np.array_equal(b,b_old):
				step=step-1
				break

			#fit again
			f=fit(p,ye,b&global_mask)
			s=np.std(p[-1][b]-f[b])
			b_old=b

		if verbose:
			print('')
			print('FITTING RESULTS:')
			print('Number of iterations requested:    ',n)
			print('Number of iterations performed:    ', step+1)
			print('Initial standard deviation:        ', s_initial)
			print('Final standard deviation:          ', s)
			print('Number of rejected points:         ',len(np.invert(b[np.invert(b)])))
			print('')
		
		return f,tmp_results,b


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
		ccf_no_edges=ccf
		ccf_no_edges[:4]=0
		ccf_no_edges[-4:]=0
		max_index=np.argmax(ccf_no_edges)
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


def galah_kern(fwhm, b):
	""" Returns a normalized 1D kernel as is used for GALAH resolution profile """
	size=2*(fwhm/2.355)**2
	size_grid = int(size) # we limit the size of kernel, so it is as small as possible (or minimal size) for faster calculations
	if size_grid<7: size_grid=7
	x= scipy.mgrid[-size_grid:size_grid+1]
	g = scipy.exp(-0.693147*np.power(abs(2*x/fwhm), b))
	return g / np.sum(g)

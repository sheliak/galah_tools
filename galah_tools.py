import pyfits
import numpy as np
import psycopg2 as mdb
import csv
import os
import fnmatch
import cPickle as pickle
from scipy import spatial
import ftputil
import getpass
import copy
from sclip.sclip import sclip
from scipy.interpolate import UnivariateSpline
import scipy
from scipy import signal
from pyflann import *

class spectrum:
	def __init__(self, name, kind='norm', extension=4, wavelength='default', linearize=True):
		
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
			path=setup.folder+str(self.date)+'/combined/'+self.name+'.fits'
		else:
			path=setup.folder+self.name+'.fits'
		
		try:#read if it exists
			hdulist = pyfits.open(path)
		except:#otherwise download
			if setup.download:
				print ' - Spectrum %s not found. Serching in downloaded spectra.' % self.name
				try:
					path=setup.folder+self.name+'.fits'
					hdulist = pyfits.open(path)
					print ' + Spectrum %s already downloaded.' % self.name
				except:
					print ' - Spectrum %s not found. Downloading from the ftp.' % self.name
					try:
						with ftputil.FTPHost('site-ftp.aao.gov.au', 'galah', getpass.getpass()) as host:
							if self.combine_method>=1:
								host.download('reductions/Iraf_5.0/%s/combined/%s.fits' % (self.date, self.name), setup.folder+self.name+'.fits')
							else:
								host.download('reductions/Iraf_5.0/%s/individual/%s.fits' % (self.date, self.name), setup.folder+self.name+'.fits')
							path=setup.folder+self.name+'.fits'
							hdulist = pyfits.open(path)
						print ' + Spectrum %s succesfully downloaded.' % self.name
					except:
						print ' + Spectrum %s failed to download.' % self.name
			else:
				print ' - Spectrum %s not found. Enable download to get it from the ftp site.' % self.name

		#set l, f, and fe
		instance={'norm':4, 'normalized':4, 'flux':0, 'fluxed':0}

		#set radial velocity if it will be needed in the future:
		#if is here because reading a spectrum is faster if read in its original velocity frame
		if (wavelength=='observer' and instance[kind]==4) or (wavelength=='object' and instance[kind]<4) or instance[kind]==4:
			con=setup.con
			if con!='':
				cur=con.cursor()
				cur.execute("select v from iraf_dr51 where name=%s" % self.name)
				try:
					self.v=float(cur.fetchone()[0])
				except TypeError:
					print ' ! Warning: no velocity in the database. Assuming v=0.'
					self.v=0.0
			else:
				self.v=float(setup.db_dict[self.name]['v'])

		try:
			self.f=hdulist[instance[kind]].data
			if instance[kind]==4:
				self.fe=hdulist[1].data
			else:
				self.fe=hdulist[instance[kind]+1].data
			crval=hdulist[instance[kind]].header['CRVAL1']
			crdel=hdulist[instance[kind]].header['CDELT1']
			self.l=np.linspace(crval, crval+crdel*len(self.f), len(self.f))
			
			if instance[kind]==4:
				#because normalized spec doesn't has its error, we use the fluxed error, but have to shift and interpolate it to normalized l:
				crval=hdulist[1].header['CRVAL1']
				crdel=hdulist[1].header['CDELT1']
				naxis=hdulist[1].header['NAXIS1']
				error_l=np.linspace(crval, crval+crdel*naxis, naxis)
				error_l=error_l*(1-self.v/299792.458)
				self.fe=np.interp(self.l,error_l,self.fe)

			
		except:
			raise RuntimeError('Cannot read spectrum. Fits extension might be missing.')
		
		#shift into correct velocity frame
		if wavelength=='default':
			pass
		elif wavelength=='observer' and instance[kind]==4:
			self.l=self.l*(1+self.v/299792.458)
		elif wavelength=='object' and instance[kind]<4:
			self.l=self.l*(1-self.v/299792.458)
		else:
			pass

		#linearize
		if linearize==True:
			self.f=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.f)
			self.fe=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.fe)
		else:
			pass
	
	def linearize(self):
		"""
		take whatever the sampling is and linearize it
		"""
		self.f=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.f)
		self.fe=np.interp(np.linspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.fe)
	
	def logarize(self):
		"""
		take whatever the sampling is and make a logaritmic sampling
		"""
		self.f=np.interp(np.logspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.f)
		self.fe=np.interp(np.logspace(self.l[0],self.l[-1],num=len(self.l)),self.l,self.fe)
	
	def shift(self, rv, linearize=True):
		"""
		shift the spectrum for radial velocity rv, given in km/s
		if linearize=True, the returned wavelength array will be linearized and flux interpolated
		"""
		l=self.l*(1-self.v/299792.458)
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

		#check if the correct resolution map is already in the memory:
		if self.ccd==1:
			if resolution_maps.map_ccd_1==False:
				resolution_maps(self.ccd)
				map_l,map_f=resolution_maps.map_ccd_1
			else:
				map_l,map_f=resolution_maps.map_ccd_1
		if self.ccd==2:
			if resolution_maps.map_ccd_2==False:
				resolution_maps(self.ccd)
				map_l,map_f=resolution_maps.map_ccd_2
			else:
				map_l,map_f=resolution_maps.map_ccd_2
		if self.ccd==3:
			if resolution_maps.map_ccd_3==False:
				resolution_maps(self.ccd)
				map_l,map_f=resolution_maps.map_ccd_3
			else:
				map_l,map_f=resolution_maps.map_ccd_3
		if self.ccd==4:
			if resolution_maps.map_ccd_4==False:
				resolution_maps(self.ccd)
				map_l,map_f=resolution_maps.map_ccd_4
			else:
				map_l,map_f=resolution_maps.map_ccd_4

		#check if wavelength calibration is linear:
		if (self.l[1]-self.l[0])==(self.l[-1]-self.l[-2]):
			linear=True
		else:
			l=self.l
			self.linearize()
			linear=False

		#extract the correct pivot number from the map:
		map_f=map_f[self.pivot-1]

		#current sampling:
		sampl=self.l[1]-self.l[0]

		#target sigma coresponds to the R=22000. We want the constant sigma, not constant R, so take the sigma that corresponds to the average R=22000
		s_target=np.ones(len(map_l))*np.average(map_l)/22000.

		#the sigma of the kernel is:
		s=np.sqrt(s_target**2-np.divide(map_l,map_f)**2)

		#fit it with the polynomial, so we have a function instead of sampled values:
		map_fit=np.poly1d(np.polyfit(map_l, s/sampl, deg=6))

		#create an array with new sampling. The first point is the same as in the spectrum:
		l_new=[map_l[0]]

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


	def knn(self,method='FLANN', K=10, d='euclidean', windows='', pickle_folder='pickled_spectra'):
		"""
		find nearest neighbours. spectra2pickle must be run first
		windows is a filename with the description of windows to use or a 1D np.ndarray
		"""

		if pspectra.names==0:
			pspectra(pickle_folder)

		spectra=pspectra.spectra
		names=pspectra.names
		l=pspectra.space

		f=np.interp(l,self.l,self.f)

		#check if window function already exists:
		if type(windows).__name__=='ndarray':#if given as an array only check if format is right
			if type(window).__name__!='ndarray': 
				raise RuntimeError('windows is not type numpy.ndarray')
			if len(window)!=len(l):
				raise RuntimeError('windows has wrong size')
			window=windows
			window_function.window=windows
			window_function.l=l
		elif len(window_function.window)==0:#if it doesnt exist create it
			window_function(l, windows)
			window=window_function.window
		else:#if it exists just use it
			window=window_function.window

		#create a mask where window==0, because there we don't have to comapre the vectors, because the difference is always 0:
		mask=np.array(window,dtype=bool)


		if method=='FLANN':
			distance={'manhattan': 'manhattan', 'euclidean': 'euclidean'}

			if flann_index.flann==False: 
				flann_index(spectra[:,mask]*window[mask], distance[d])
			
			ind,dist=flann_index.flann.nn_index(f[mask]*window[mask],K,checks=flann_index.index['checks'])

			if distance[d]=='euclidean':
				return names[ind], np.sqrt(dist)
			else:
				return names[ind], dist

		if method=='KDTree':
			distance={'manhattan': 1, 'euclidean': 2}

			if kdtree_index.index==False:
				kdtree_index(spectra[:,mask]*window[mask])

			dist,ind=kdtree_index.index.query(f[mask]*window[mask],K,p=distance[d])

			return names[ind], dist

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

class pspectra:
	spectra=0
	names=0
	space=[]

	def __init__(self,pickle_folder):
		if pickle_folder[-1]=='/': pickle_folder=pickle_folder[:-1]

		try:
			pspectra.space=pickle.load(open('%s/space.pickle' % (pickle_folder), 'rb'))
			#f=np.interp(space,self.l,self.f)
			#l=space
		except:
			raise RuntimeError('Pickle spectra first if you want to use nearest neighbour search. No wavelength samling is saved.')

		blocks=[]
		for path, dirs, files in os.walk(os.path.abspath(pickle_folder)):
			for filename in fnmatch.filter(files, 'b*.pickle'):
				blocks.append(filename)

		if len(blocks)==0:
			raise RuntimeError('Pickle spectra first if you want to use nearest neighbour search. No pickled spectra saved.')

		spectra_frompickle=[]
		for i in blocks:
			spectra_frompickle=spectra_frompickle+pickle.load(open('%s/%s' % (pickle_folder,i), 'rb'))

		names=[]
		spectra=[]

		for i in spectra_frompickle:
			names.append(i[0])
			spectra.append(i[1])

		pspectra.spectra=np.array(spectra)
		pspectra.names=np.array(names)

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

class flann_index:
	flann=False
	index=0

	def __init__(self, spectra, d):
		set_distance_type(d)
		flann_index.flann=FLANN()
		flann_index.index=flann_index.flann.build_index(spectra, algorithm='autotuned', target_precision=0.9)

	def clear(self):
		flann_index.flann=False
		flann_index.index=0

class kdtree_index:
	index=False

	def __init__(self, spectra):
		kdtree_index.index=spatial.cKDTree(spectra)

	def clear(self):
		kdtree_index.index=False

class resolution_maps:
	map_ccd_1=False
	map_ccd_2=False
	map_ccd_3=False
	map_ccd_4=False

	def __init__(self, ccd):
		hdulist = pyfits.open('resolution_maps/ccd%s_piv.fits' % ccd)
		res=hdulist[0].data
		crval=hdulist[0].header['CRVAL1']
		crdel=hdulist[0].header['CDELT1']
		naxis=hdulist[0].header['NAXIS1']
		l=np.linspace(crval, crval+crdel*naxis, naxis)
		hdulist.close()

		if ccd==1:
			resolution_maps.map_ccd_1=(l,res)
		if ccd==2:
			resolution_maps.map_ccd_2=(l,res)
		if ccd==3:
			resolution_maps.map_ccd_3=(l,res)
		if ccd==3:
			resolution_maps.map_ccd_3=(l,res)


class setup:
	folder=''
	folder_is_root=False
	con=''
	csv=''
	db_dict={}
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

			if key=='con':
				setup.con=kwargs['con']

			if key=='download':
				setup.download=kwargs['download']

			if key=='csv':
				setup.csv=kwargs['csv']
				reader = csv.DictReader(open(setup.csv))
				for row in reader:
					key = row.pop('name')
					if key in setup.db_dict:
						pass
					setup.db_dict[key] = row

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


def spectra2pickle(ccd, space=[], limit=999999999999, pickle_folder='pickled_spectra'):
	"""
	Translate all spectra into blocks of pickle objects. Blocks are used, because only 10000 spectra can be in one object. 
	After the pickle blocks are read they can be combined into a single table again.
	This is used by nearest neighbour search, because spectra must be read fast.
	"""

	if pickle_folder[-1]=='/': pickle_folder=pickle_folder[:-1]

	#using pickle to store data can be dangerous. Display warning
	print ' ! Warning: using pickle to store data can be dangerous. Use at your own risk!'

	#check if defined space is valid for the given ccd:
	if ccd==1:
		if space==[]: space=np.linspace(4713,4903,4100)
		if not any(4713<i<4903 for i in space): raise RuntimeError('Space is out of boundary for ccd 1')
	if ccd==2:
		if space==[]: space=np.linspace(5600,5875,4100)
		if not any(5600<i<5875 for i in space): raise RuntimeError('Space is out of boundary for ccd 2')
	if ccd==3:
		if space==[]: space=np.linspace(6477,6740,4100)
		if not any(6477<i<6740 for i in space): raise RuntimeError('Space is out of boundary for ccd 3')
	if ccd==4:
		if space==[]: space=np.linspace(7584,7887,4100)
		if not any(7584<i<7887 for i in space): raise RuntimeError('Space is out of boundary for ccd 4')

	#create a folder to store pickled spectra:
	if not os.path.exists('pickled_spectra'):
		os.makedirs('pickled_spectra')

	#create a list of spectra:
	spectra=[]
	for path, dirs, files in os.walk(os.path.abspath(setup.folder)):
		for filename in fnmatch.filter(files, '*%s.fits' % (ccd)):
			spectra.append(filename[:-5])
	
	#read and interpolate spectra one by one:
	block=[]
	block_num=0
	n=0
	nn=0
	for i in spectra:
		if n==9999:#only 10000 spectra can be saved in one file
			pickle.dump(block,open('%s/b%s.pickle' % (pickle_folder,block_num), 'wb'))
			block_num+=1
			n=0
			block=[]
		try:#if there is no normalized spectrum skip it
			s=read(i, kind='norm', wavelength='default').interpolate(space)
			block.append([i,s.f])
			n+=1
		except RuntimeError:
			pass
		nn+=1

		if nn>=limit: break# stop reading spectra if limit is reached (to spare the memory or time when testing)

		if nn%10 == 0: print nn, '/', len(spectra), "pickled."

	
	pickle.dump(block,open('%s/b%s.pickle' % (pickle_folder,block_num+1), 'wb'))
	pickle.dump(space,open('%s/space.pickle' % (pickle_folder), 'wb'))

	return space

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

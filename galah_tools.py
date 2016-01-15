import pyfits
import numpy as np
import psycopg2 as mdb
import csv
import os
import fnmatch
import cPickle as pickle
from pyflann import *
from scipy import spatial
import ftputil
import getpass

class spectrum:
	def __init__(self, name, kind='norm', extension=4, wavelength='default', linearize=True):
		
		#set atributes
		self.name = name
		self.ccd=int(self.name[-1])
		self.date=int(self.name[:6])
		self.run=int(self.name[6:10])
		self.combine_method=int(self.name[10:12])
		self.pivot=int(self.name[12:15])

		#read spectrum
		self.l=-1
		self.f=-1

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
		
		instance={'norm':4, 'normalized':4, 'flux':0, 'fluxed':0}

		try:
			self.f=hdulist[instance[kind]].data
			crval=hdulist[instance[kind]].header['CRVAL1']
			crdel=hdulist[instance[kind]].header['CDELT1']
			self.l=np.linspace(crval, crval+crdel*len(self.f), len(self.f))
		except:
			raise RuntimeError('Cannot read spectrum. Fits extension might be missing.')

		#set velocity frame
		con=setup.con
		if con!='':
			cur=con.cursor()
			cur.execute("select v from iraf_dr50 where name=%s" % self.name)
			try:
				self.v=float(cur.fetchone()[0])
			except TypeError:
				print ' ! Warning: no velocity in the database. Assuming v=0.'
				self.v=0.0
		else:
			self.v=float(setup.db_dict[self.name]['v'])

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
		else:
			pass

	def shift(self, rv, linearize=True):
		"""
		shift the spectrum for radial velocity rv, given in km/s
		if linearize=True, the returned wavelength array will be linearized and flux interpolated
		"""
		l=self.l*(1-self.v/299792.458)
		if linearize==True:
			self.f=np.interp(np.linspace(l[0],l[-1],num=len(l)),l,self.f)
			self.l=np.linspace(l[0],l[-1],num=len(l))
		else:
			self.l=l

		return self

	def add_noise(self,noise):
		pass

	def interpolate(self,space):
		"""
		interpolate the spectrum to a wavelength space defined in space. 
		"""
		space=np.array(space)
		self.f=np.interp(space,self.l,self.f)
		self.l=space
		return self

	def normalize(self):
		pass

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


class flann_index:
	flann=False
	index=0

	def __init__(self, spectra, d):
		set_distance_type(d)
		flann_index.flann=FLANN()
		flann_index.index=flann_index.flann.build_index(spectra, algorithm='autotuned', target_precision=0.7)

	def clear(self):
		flann_index.flann=False
		flann_index.index=0

class kdtree_index:
	index=False

	def __init__(self, spectra):
		kdtree_index.index=spatial.cKDTree(spectra)

	def clear(self):
		kdtree_index.index=False


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


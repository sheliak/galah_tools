#This is a script that shows the nearest neighbour search possibilities.

#Import galah tools:
import galah_tools as gtools

#Import a module for postgresql, if you want to use it (see next step):
import psycopg2 as mdb

#Chose one way (only one) of setting the working environment. This software comes without a database, because it is public. The database must be downloaded from a proper GALAH webpage/ ftp site.

#1. Here we set the folder where data is in standard folders and a connection to the sql database
gtools.setup(root_folder='./test/iraf_dr51', con=mdb.connect("dbname=hermes_master user=janez"), download=False)

#2. Here we set a folder where all the spectra are collected without any hierarchy. We set the csv file as a source of database. 
#gtools.setup(folder='/home/janez/fax/galah/tools/test/random_spectra', csv='test/db/iraf_dr50.csv')


#Read a spectrum. Several options are available. By default a normalized spectrum will be read in whatever samping it is given in the fits file. The spectrum will be called s1:
s1=gtools.read('1402070012010053')

#We can force the spectrum to be linearized and in the star's velocity frame (spectral lines have laboratory wavelengths).
s2=gtools.read('1402070012010053', linearize=True, wavelength='object')

#If spectrum is not in the given folders, it can be downloaded if download option is enabled in setup. 
s3=gtools.read('1402070012010063')

#Both spectra have properties. .l and .f return the wavelength sampling and fluxes for this sampling as arrays.
print s1, s1.l, s1.f
print s2, s2.l, s2.f

#you can print some properties for the spectra
print s1, s1.name, s1.ccd, s1.date ,s1.run, s1.combine_method, s1.pivot

#we can shift the spectrum for 20 km/s
#s1= s1.shift(20)

#pickle all the spectra if you want to use nearest neighbour search. Pickling the spectra first is mandatory, otherwise reading them again and again will be too slow.
#here we do it for ccd 3
#returned array l is the default sampling that is hard coded in the code. You can define an arbitrary sampling with space keyword (must be a 1D array).
#pickling a large number of spectra is time consuming. Don't do it each time. The old pickled spectra will be saved in a folder.
#A folder where pickled spectra will be saved can be given with pickle_folder keyword. Default is ./pickled_spectra
l=gtools.spectra2pickle(3)

#We will use some windows when we do the nearest neighbour search. Let us plot the windows first:
w=gtools.window_function(l, 'test/windows.txt')
w.plot_window()

#Here we clear the saved windows. This step is only necessary if new windows will be defined in the next step, so it can be skipped in this example.
w.clear()

#find 2 nearest neighbours for the spectrum s1. We use a file with defined windows.
#A folder where pickled spectra are saved can be given with pickle_folder keyword. Default is ./pickled_spectra
print s1.knn(K=2, windows='test/windows.txt')

#find 2 nearest neighbours for the spectrum s2. The same windows will be used as before. Software saves as many variables as possible, like windows, so they don't have to be recalculated in each step.
print s1.knn(K=2, method='KDTree')

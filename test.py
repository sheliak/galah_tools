#This is a script that shows the nearest neighbour search possibilities.

#Import galah tools:
import galah_tools as gtools
from matplotlib import *
from pylab import *


#Chose one way (only one) of setting the working environment. This software comes without a database, because it is public. The database must be downloaded from a proper GALAH webpage/ ftp site.

#1. Here we set the folder where data is in standard folders and a connection to the sql database
gtools.setup(root_folder='../obdelava/reductions/results')

#2. Here we set a folder where all the spectra are collected without any hierarchy. We set the csv file as a source of database. 
#gtools.setup(folder='/home/janez/fax/galah/tools/test/random_spectra', csv='test/db/iraf_dr50.csv')


#Read a spectrum. Several options are available. By default a normalized spectrum will be read in whatever samping it is given in the fits file. The spectrum will be called s1:
s1=gtools.read('1707110015010062')

#We can force the spectrum to be linearized and in the star's velocity frame (spectral lines have laboratory wavelengths).
s2=gtools.read('1707110015010063', linearize=True, wavelength='object')

#If spectrum is not in the given folders, it can be downloaded if download option is enabled in setup. 
#s3=gtools.read('1807110017000062')

#Both spectra have properties. .l and .f return the wavelength sampling and fluxes for this sampling as arrays.
print s1, s1.l, s1.f
print s2, s2.l, s2.f

#you can print some properties for the spectra
print s1, s1.name, s1.ccd, s1.date ,s1.run, s1.combine_method, s1.pivot

plot(s1.l,s1.f,'k-')

s1.equalize_resolution()
#s1.convolve(0.22)
#s1.res_degradation(22000,15000)

plot(s1.l,s1.f,'r-')

show()

#we can shift the spectrum for 20 km/s
#s1= s1.shift(20)

#We will use some windows when we do the nearest neighbour search. Let us plot the windows first:
#w=gtools.window_function(l, 'test/windows.txt')
#w.plot_window()

#Here we clear the saved windows. This step is only necessary if new windows will be defined in the next step, so it can be skipped in this example.
#w.clear()

#find 2 nearest neighbours for the spectrum s1. We use a file with defined windows.
#A folder where pickled spectra are saved can be given with pickle_folder keyword. Default is ./pickled_spectra
#print s1.knn(K=2, windows='test/windows.txt')

#find 2 nearest neighbours for the spectrum s2. The same windows will be used as before. Software saves as many variables as possible, like windows, so they don't have to be recalculated in each step.
#print s1.knn(K=2, method='KDTree')

#GALAH TOOLS PACKAGE
Galah tools package is a python module enabling an easy acces to GALAH spectra and some products. It also includes several useful tools for manipulating and analysing the spectra.

##Download

This software can be downloaded by cloning the git repository:

```{r, engine='bash'}
git clone https://github.com/sheliak/galah_tools.git
```

A folder `galah_tools` will be created

##Instalation

A setup.py file is provided. The module can be installed in the usual way by running

```{r, engine='bash'}
cd galah_tools
sudo python setup.py install 
```
This will make the `galah_tools` module available system-wide.

##Test

A short test scrip is provided to test the basic funcionality and the search for nearest neighbours. See the contest of the `test.py` and run it with

 ```{r, engine='bash'}
python test.py
```

See the contest of `test/windows.txt` to learn how to set up windows and ranges used in nearest neighbour search.


#Usage instructions

##Setup of the environment

Before you start exploiting the `galah_tools` module, you have to set the working environment. This way the module will know where to look for the spectra and where to look for the data tables. 

This software comes with no data tables and with spectra of only two random stars, because it is publicly available, unlike the spectra or tables. You will need the acces to GALAH data to use this software.

The environment is set up by calling the setup function. At least two parameters must be given: folder where the spectra are stored and connection to the database.

The following example set-ups the environment where the spectra are stored in a single folder with no hirearchy and the csv dump is used as a database.

 ```python
import galah_tools as gtools

gtools.setup(folder='./test/random_spectra', csv='test/db/iraf_dr50.csv')
```

The second option is to use the same folder structure as on the ftp site for the spectra storage and an sql database instead of csv table. This will also be faster. Any combination of these options can be used.

 ```python
import galah_tools as gtools
import psycopg2 as mdb #this is the preffered module to use with the postgre sql database

gtools.setup(root_folder='./test/iraf_dr50', con=mdb.connect("dbname=hermes_master user=janez"))
```
Notice that instead of `folder` parameter we used `root_folder`. 

If you want to use spectra that are not saved in your local system, you can give `galah_tools` a permission to download them on the go. They will be downloaded from the ftp site and you will be prompted for the password. 

```python
import galah_tools as gtools

gtools.setup(folder='./', csv='test/db/iraf_dr50.csv', download=True)
```

##Reading the spectra

Spectrum can be read by calling its name (a string or an int can be used)

```python
s=gtools.read('1402070012010053')
```

Several arguments can be given to control the properties of the spectrum:
* `kind={'norm{alized}'} {'flux{ed}'}} {'nosky'}` to retrieve the normalized or fluxed spectrum or the spectrum before the sky subtraction.
* `extension={0}{1}{2}{3}{4}` to retrieve a spectrum from the specific fits extension. 0=fluxed, 2=spectrum before the sky subtraction and 4=normalized.
* `linearize={True} {False}` If you want to force the spectrum to have a linearized sampling set to True. Otherwise you will get whatever is written in the fits file.
* `log={True} {False}` If you want to use the log spacing in the wavelength set this to `True`. **Not yet implemented. Do not use.**
* `wavelength={'default'} {'observer'}, {'object'}, {'bary'}` This controls in what velocity space you want the spectrum. `'default'` gives whatever is in the fits file. `'observer'` returns a spectrum with the wavelengths as measured from the arc lamp, uncorrected for the barucentric velocity. `'object'` gives the wavelengths corrected for the barycentric velocity and the radial velocity of the star (as measured by GUESS). `'bary'` gives the wavelengths corrected for the barycentric velocity. **Please, use only `'object'` for now. Other options have not been tested yet or are in the development**

When spectrum is downloaded from the GALAH ftp site it is saved in the `folder` or `root_folder`, whichever is defined.

##Spectrum class

The class has the following attributes:
*`s.l` is the wavelength array
*`s.f` is the flux array
*`s.name` is the name of the spectrum (the same one you use to read the spectrum)
*`s.ccd` is the ccd number for this spectrum
*`s.date` is the yymmdd date mark
*`s.run` is the run number
*`s.combine_method` is the combine metod (integer)
*`s.pivot` is the pivot number
*`s.v` is the radial velocity (as measured by GUESS)

##Modyfing spectra

Spectra can be quickly modified (shifted for a specified v, normalized, etc.). A universal synatx is:
```python
s.modifier(args)
```
This will rewrite the attributes of the opened spectrum. If you want to save the modified spectrum into a different variable, you have to hard copy it. 

###Normalization
**TO DO**

###Radial velocity shift
```python
s.shift(rv, linearize=True)
```
Shift for `rv`, given in km/s. If `linearize=True`, linearize the wavelength scale. If False it will remain as is after the rv shift is applied, regardless the `linearize` parameter when opening the spectra in the first place.

###Add noise
```python
s.add_noise(snr, target_snr, skip=True)
```
Change the SNR of the spectrum from `snr` to `target_snr`. If `target_snr>=snr`, `skip` is checked. If True, nothing will happen and the SNR will remain as is. If set to False an error will pop up.

**To do: read the original snr from the database or from the error spectrum. Give option to decide which snr to use.**

###Interpolate
```python
s.interpolate(space)
```
Interpolate the spectrum into a given `space`. `space` can be any 1D list or array filled with numeric values. `s.l` will be changed into space and `s.f` will be lineary interpolated.

##Other methods

##Finding nearest neighbours

There are currently two methods available for finding nearest neighbours; scipy's KNN and FLANN. 

###Pickling spectra

It is mandatory to pickle the spectra before the nearest neighbour search is called. Pickling spectra means that they will be read into the memory and the correct part of the memory will be saved to the hard drive in a pickle format. Reading the spectra initially takes around an hour for all the GALAH spectra and reading back the pickled spectra takes minutes. So from there on, you will only have to wait minutes every time you run the nearest neighbour search instead of an hour. If you run several nearest neighbour searches in one script, the pickled spectra will be read only ones and saved in the memory until the script is terminated.

Spectra can be pickled by running
```python
l=gtools.spectra2pickle(ccd, space, limit, pickle_folder)
```
Only the first argument is mandatory, telling spectra for which arm to pickle. `space` can be a wavelength sampling into which the spectra will be interpolated, `limit` is a maximum number of spectra to pickle, usable for testing pourposes, and `pickle_folder` is the name of the folder where pickled spectra will be saved. By default this is `./pickled_spectra`. The wavelength space `l` is returned if pickling is successful.

10000 spectra will be saved in the same pickle object. If there are more spectra to pickle, they will be saved into several files. In addition, the wavelength space is saved too.

**Pickled objects should never be sent over the internet or exchanged with untrustworthy people, because there is no safety check when they are read back into the memory!**

###Defining windows and ranges

When looking for nearest neighbours, you want to give different weights to different parts of the spectrum. We provide a simple tool for defining the weights in a text file. See `test/windows.txt` for the description and examples.

Weights can also be given as a 1D numpy array of the same length as the spectra are. 

If a weight at a certain pixel is zero, the pixel will be filtered out and the nearest neighbour search will be able to work with fewer dimensions, so it will be faster. Consider this when setting weights to values close to zero.

Windows can be retrieved from the file and, for example, plottted:
```python
w=gtools.window_function(l, 'test/windows.txt')
w.plot_window()
```

###Running the nearest neighbour search
```python
names, distances=s.knn(K, windows, method, d, pickle_folder)
```
Ony the first argument is mandatory.
* `K` is the number of neighbours you want (must be smaller or equal to number of spectra)
* `windows` is a filename where windows and ranges are defined
* `method={KDTRee}{FLANN}` tells which method to use
* `d={manhattan}{euclidean}` tells what metric to use
* `pickle_folder` is the name of the folder where pickled spectra are saved. Default is `./pickled_spectra`
* `names` and `distances` are arrays with the nearest neighbours arranged by distance. First one gives names of the nearest spectra and the second one gives distances.


#Licence

Copyright (C) 2015  Janez Kos

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

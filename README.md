#GALAH TOOLS PACKAGE
Galah tools package is a python module enabling an easy acces to GALAH spectra and some products. It also includes several useful tools for manipulating and analysing the spectra.

##Instalation

A setup.py file is provided. The module can be installed in the usual way by running

```{r, engine='bash'}
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
* `linearize={True}{False}` If you want to force the spectrum to have a linearized sampling set to True. Otherwise you will get whatever is written in the fits file.
* `wavelength={'default'} {'observer'}, {'object'}, {'bary'}` This controls in what velocity space you want the spectrum. `'default'` gives whatever is in the fits file. `'observer'` returns a spectrum with the wavelengths as measured from the arc lamp, uncorrected for the barucentric velocity. `'object'` gives the wavelengths corrected for the barycentric velocity and the radial velocity of the star (as measured bu GUESS). `'bary'` gives the wavelengths corrected for the barycentric velocity. **Please, use only `'object'` for now. Other options have not been tested yet or are in the development**


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

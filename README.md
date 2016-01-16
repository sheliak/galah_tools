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

 ```{r, engine='python'}
python test.py 
```

See the contest of `test/windows.txt` to learn how to set up windows and ranges used in nearest neighbour search.


#Usage instructions

##Setup of the environment

Before you start exploiting the `galah_tools` module, you have to set the working environment. This way the module will know where to look for the spectra and where to look for the data tables. 

This software comes with no data tables and with spectra of only two random stars, because it is publicly available, unlike the spectra or tables. You will need the acces to GALAH data to use this software.

The environment is set up by calling the setup function. At least two parameters must be given: folder where the spectra are stored and connection to the database.

###Spectra storage

Spectra can be all stored in a single folder with no hirearchy

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

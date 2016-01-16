#GALAH TOOLS PACKAGE
Galah tools package is a python module enabling an easy acces to GALAH spectra and some products. It also includes several useful tools for manipulating and analysing the spectra.

##Instalation

A setup.py file is provided. It can be installed in the usual way by running

```{r, engine='bash'}
sudo python setup.py install 
```
This will make the `galah_tools` module available system-wide.

##Test

* A short test scrip is provided to test the basic funcionality and the search for nearest neighbours. See the contest of the `test.py` and run it with

 ```{r, engine='python'}
python test.py 
```

* See the contest of `test/windows.txt` to learn how to set up windows and ranges used in nearest neighbour search.


#Usage instructions

##Setup of the environment

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

Zebrafish scale image processing, Erk activity and tissue flows quantification sample code from:
Alessandro De Simone, Maya N. Evanitsky, Luke Hayden, Ben D. Cox, Julia Wang, Valerie A. Tornini, Jianhong Ou, Anna Chao, Kenneth D. Poss, Stefano Di Talia
Control of osteoblast regeneration by a train of Erk activity waves 

Copyright (C) 2020  Alessandro De Simone

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


This code will perform image processing to quantify Erk activity and
tissue flows from images of regenerating zebrafish scales.

The code requires MATLAB (MathWorks).
The code runs in MATLAB_R2016b. The code runs on a 
iMac (27-inch, Late 2013) macOS Sierra
Processor 3.5 GHz Intel Core i7
Memory 32 GB 1600 MHz DDR3

INSTALLATION AND USAGE
The user can run the code by changing 
MATLAB's present working directory to 'DeSimoneErkwaves2020-master'.
MS Windows users have to adapt paths to MS Windows sintax.
 
Data files are too large to include within this Github repo and must be
downloaded separately from:
https://drive.google.com/file/d/117tJhklJNEQJ3ND4D76V_xl22kcwfzZ4/view?usp=sharing
The downloaded folder must be unzipped and tne 'data' folder must be
placed in 'DeSimoneErkwaves2020-master'.

MATLAB Toolboxes required are:
- Image Processing Toolbox (version 9.5)
- Statistics and Machine Learning Toolbox (version 11.0)
- Financial Toolbox (version 5.8)
- Curve Fitting Toolbox (3.5.4)

The code requires external MATLAB functions:
loadtiff.m, saveastiff.m, weightedcov.m, xml2struct.m, nanconv.m,
brewermap.m
Those functions be downloaded from MathWorks. We attach them to the code together with
their licenses.

The main image processing routine is the script
CREATE_H2A_ERKKTR_EXAMPLE.M
The script can run altogether or each section sequentially.
In some instances, user input would be required, but
saved parameters are used for demo purposes. Users wishing to 
use custom parameters can change USESAVEDPARAMS to false
in the appropriate sections.

Sample outputs for each steps are provided in the data folder.

Scale data are saved in objects files in the "objects" folder (myScale
data struct). Stacks are saved as tiff files.

Variables that an user could change are briefly
described. 

DATA SAMPLE

A sample dataset is provided is provided (see Installation). Erk activity and tissue flow maps can
be generated using the code. 

Time-points include two channels:
ch1 -> Erk sensor (osx:ERKKtr-mCerulean)
ch2 -> osteoblast nuclear marker (osx:H2A-mCherry)

Expected output for each step is provided. 
16 Gb RAM is required
Data sample runtime: ~3h

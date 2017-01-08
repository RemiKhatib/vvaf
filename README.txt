#This program has been written by Rémi Khatib.
#It allows you to calculate several kind of
#velocity-velocity correlation functions.
#The goal is to simulate vibrational spectra
#(IR, Rman, SFG)


#========
#Makefile
#========
#Compilation
> make

#Clean the obj directory
> make clean

#Clean everything created during the compilation
> make very_clean

#Create a TAGS file for emacs (etag)
> make tag

#You may change the compilation options by changing DEBUG in the
#makefile
DEBUG='opt'     #Will just perform an O2 optimization
DEBUG='warning' #Will use some common warnings
DEBUG='fprof'   #Same than warning but with a fprofile

#========
#Tutorial
#========
The methodology is described in __benchmarks. Please,
follow the numbers.
In each directory, you will find an input file where
the keywords will be explained and a command.sh file
to lauch directly the calculation of the correlation
functions

To facilitate the reading of the input files, you
may use the plugins for emacs or vi
https://www.cp2k.org/tools


#=========
#Equations
#=========
If you want to read few equations about first order expansion
, or how to go from the bond framework to the molecular
framework, please refer to:
./__1st_expansion/equation.pdf


#====
#Bugs
#====
There are some minor bugs (the layering loop until the infiny
for some parameters associated with the slabs, ...). But
the main features should be good.
Please contact me if there is any problem.

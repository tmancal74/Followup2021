This file contains information about simulation scripts for the manuscript:

 Veronica R. Policht, Cameron Spitzfaden, Jennifer P. Ogilvie and Tomáš Mančal
 ..., submitted 2021

The following files should be present in the current directory:

README.txt  (this file)
script_Followup2021.py
script_Followup2021.yaml
Makefile
make.bat
scr/aux_figures.py
scr/aux_movies.py
scr/parula_colormap.dat
scr/probe_mpi.py

Quick reference:

SOFTWARE INSTALLATION
---------------------
The scripts require Python 3 and the Quantarhei package installed. Open the file
"script_Folloup2021.py" as a text file, and consult it for the installation
instructions of Python 3 and the Quantarhei package.

RUNNING THE SCRIPTS
-------------------
The simplest way to run the script is to configure the session in the first few
lines of the file called Makefile (this applies to both Windows and Linux/Unix/Mac
systems, although Makefile is used very differently on different systems).
Set the options depending on your software/hardware and use the following command:

make run

This again applies to all systems. On Windows Makefile is just a configuration file
for a custom script, while on Linux/Unix/Mac it is used by a unix specific program
called make.

On Windows you can alternatively click on the runme.bat file to start the scripts in
a new window. Path to you Python interpreter has to be set on your system for this
feature to work.

More information about more control over running the scripts can be found in
script_Followup2021.py file.

CONFIGURATION AND PARAMETERS OF THE SIMULATIONS
-----------------------------------------------

Consult the file "script_Followup2021.yaml" for details on the simulation
parameters and options.



All files in this directory contain descriptions of their meaning and usage. Just
open them as text files.

Last change: March 12, 2021
Author: Tomas Mancal
Email: tmancal74@gmail.com

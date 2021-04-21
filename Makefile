################################################################################
#
#  GLOBAL SIMULATION SETTINGS
#  Here you can set behavior of the simulation
#  (Details of Makefile usage below the settings)
#
################################################################################

# Number of processes to start (if > 1, mpi4py Python package and MPI have
# to be installed for the simulation to run in parallel)
NUMBER_OF_PROCESSES=1

# run in the background
BACKGROUND=False

# filename to save output to (if BACKGROUND=True the default is output.log);
LOG_FILE=

# change this to your python interpreter
PYTHON= python

################################################################################
#
#  Makefile for the simulaiton scripts of the manuscript:
#
#  Veronica R. Policht,  Cameron Spitzfaden, Jennifer P. Ogilvie
#  and Tomáš Mančal
#	 ..., submitted 2021
#
#  This file should work on most Linux/Unix systems (including macOS)
#
#  Type the following on the command line
#  (">" represents the command line prompt)
#EXCITON_SCRIPT}
#  > make help
#
#  to see the list of available tasks. The tasks can be configured in the
#  configuration section below. Unless you know precisely what you are doing,
#  do not edit anything except the configuration section
#
################################################################################

#
#   DO NOT EDIT BELOW THIS LINE
#
SCRPTNAME=script_Followup2021
SCRDIR=scr
MOVIES_SCRIP=${SCRDIR}/aux_movies.py
FIGURES_SCRIPT=${SCRDIR}/aux_figures.py
PATHWAYS_SCRIPT=${SCRDIR}/aux_pathways.py
EXCITON_SCRIPT=${SCRDIR}/aux_excitons.py

# PBS
JOBNAME=follow


# set PARALLEL depending on the number of required processes
ifeq ($(shell test ${NUMBER_OF_PROCESSES} -gt 1; echo $$?),0)
PARALLEL=True
else
PARALLEL=False
endif

#probe MPI presence if it is required
MPI_REPORT= MPI presence not tested
ifeq (${PARALLEL},True)
ifeq ($(shell ${PYTHON} ${SCRDIR}/probe_mpi.py; echo $$?),0)
MPI_REPORT= MPI probed with success!
PARALLEL=True
else
MPI_REPORT= MPI not found \(mpi4py package or MPI implementation is missing\)
PARALLEL=False
endif
endif

# set SAVE_OUTPUT to True if a LOG_FILE is set
ifeq (${LOG_FILE},)
SAVE_OUTPUT=False
else
SAVE_OUTPUT=True
endif

PARALLELOPT=
ifeq (${PARALLEL},True)
PARALLELOPT= -p -n ${NUMBER_OF_PROCESSES}
endif

LOGGING=
ifeq (${SAVE_OUTPUT},True)
LOGGING= > ${LOG_FILE} 2>&1
endif

AMPRS=
ifeq (${BACKGROUND},True)
AMPRS=&
ifeq (${LOGGING},)
LOGGING= > output.log 2>&1
endif
endif

PIPE= ${LOGGING} ${AMPRS}

#
################################################################################
#
#   AVAILABLE TASKS
#
################################################################################
#

# default task
all: help
	@echo
	@echo Current settings:
	@echo -----------------
	@echo
	@echo NUMBER_OF_PROCESSES=${NUMBER_OF_PROCESSES}
	@echo ${MPI_REPORT}
	@echo BACKGROUND=${BACKGROUND}
	@echo LOG_FILE=${LOG_FILE}
	@echo PYTHON=${PYTHON}
	@echo
	@echo Will run in parallel: \(PARALLEL=\) ${PARALLEL}
	@echo Output will be saved: \(SAVE_OUTPUT=\) ${SAVE_OUTPUT}
	@echo

# help message
help:
	@echo
	@echo "Simulation Makefile          "
	@echo "==================="
	@echo
	@echo "To configure the session, edit the switches at the start"
	@echo "of this Makefile. You can choose between serial and parallel"
	@echo "simulation and adjust logging and whether to run in the "
	@echo "background or not "
	@echo
	@echo "Available tasks: "
	@echo "---------------- "
	@echo
	@echo "> make help "
	@echo
	@echo "    Prints this message "
	@echo
	@echo "> make run "
	@echo
	@echo "    Runs the simulations "
	@echo
	@echo "> make figures DIR=results_directory"
	@echo
	@echo "    Produces 2D omega_2 map figures from the simulation "
	@echo "    which ran in the simulation mode 'single' or 'disorder'"
	@echo "    (see configureation yaml file). results_directory is "
	@echo "    the directory containing results of Quantarhei simulation."
	@echo
	@echo "> make movies DIR=results_directory"
	@echo
	@echo "    Produces movies for the energy gap scan from the simulation "
	@echo "    which ran in the simulation mode 'scan'"
	@echo "    (see configureation yaml file). results_directory is"
	@echo "    the directory containing results of Quantarhei simulation."
	@echo
	@echo "> make pathways FILE=pws_file N=n"
	@echo
	@echo "    Prints out Feynman diagrams of rephasing pathways contributing to "
	@echo "    a spectral region specified within the script file "
	@echo "    'scr/aux_pathways.py'. N specifies the maximum number "
	@echo "    of printed pathways."
	@echo
	@echo "> make clean "
	@echo
	@echo "    Deletes the output of the simulations "
	@echo
	@echo "> make del "
	@echo
	@echo "    Deletes media files created by scripts "
	@echo

# delete results from all previous runs
clean:
	rm -rf sim* log output.log *.tar job.sh test.log

# delete media produced by auxiliary scripts
del: clean
	rm -rf *.png *.mov *.mp4 *.dat *.bak 

# delete everything
purge: clean del

# run a simulation
run:
	(time qrhei run ${PARALLELOPT} ${SCRPTNAME}.yaml) ${PIPE}


##### PBS specific tasks

job.sh: job.temp
	 @awk '{gsub("__WDIR__", wdir); gsub("__HOME__", home); gsub("__JOB_NAME__", jname); gsub("__NCPUS__", ncpus)}1' home=${HOME} wdir=`pwd` jname=${JOBNAME} ncpus=${NUMBER_OF_PROCESSES} job.temp > job.sh


qsub: job.sh
	qsub -q global@elixir-pbs.elixir-czech.cz job.sh > jobId.txt

check:
	@qstat -f `cat jobId.txt` | grep job_state
	@qstat -f `cat jobId.txt` | grep resources_used

#####

copy:
	mkdir ../${DIR}
	cp Makefile ../${DIR}
	cp script_Followup2021.py ../${DIR}
	cp script_Followup2021.yaml ../${DIR}
	cp job.temp ../${DIR}
	cp -r scr/ ../${DIR}
ifneq ("$(wildcard run.sh)","")
	cp run.sh ../${DIR}
endif

# make figures from raw data (single realization or average)
figures:
	${PYTHON} ${FIGURES_SCRIPT} ${DIR}

# make movies from raw data of enegy gap scan
movies:
	${PYTHON} ${MOVIES_SCRIP} ${DIR} ${NUMBER_OF_PROCESSES}

pathways:
	${PYTHON} ${PATHWAYS_SCRIPT} ${FILE} ${N}

excitons:
	${PYTHON} ${EXCITON_SCRIPT} ${DIR}

# creates a tar ball with all the files required to run simulations
# (Unix/Linux/macOS only feature)
pack:
	tar cf ${SCRPTNAME}.tar ${SCRDIR} make.bat Makefile
	tar rf ${SCRPTNAME}.tar README.txt runme.bat ${SCRPTNAME}.*

back:
ifneq ("$(wildcard script_Followup2021.yaml)","")
	mv script_Followup2021.yaml script_Followup2021.yaml.bak
endif

#
# Input files for example runs
#
set_example_single: back
	cp templates/script_Followup2021_example_single.yaml ./script_Followup2021.yaml
	
set_example_scan: back
	cp templates/script_Followup2021_example_scan.yaml ./script_Followup2021.yaml
	


################################################################################
# EOF

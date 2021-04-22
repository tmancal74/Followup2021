@echo off
rem #######################################################################
rem
rem     This script models the behaviour of the Unix/Linux utility
rem     called make. You do not need to install make on your Windows
rem     system to run the simulations in this directory.
rem
rem     This file uses "Makefile" as its configuration file. If you want to
rem     change the behaviour of this script, edit the opening few lines
rem     of Makefile. Makefile resides in the same directory as this script.
rem
rem #######################################################################

for /f "tokens=1,2 delims==" %%a in (Makefile) do (
if %%a==NUMBER_OF_PROCESSES set NUMBER_OF_PROCESSES=%%b
if %%a==LOG_FILE set LOG_FILE=%%b
if %%a==BACKGROUND set BACKGROUND=%%b
if %%a==PYTHON set PYTHON=%%b
)

rem ######################################################################
rem
rem     Settings (do not edit anything below this line)
rem
rem ######################################################################

set SCRPTNAME=script_Followup2021
set SCRDIR=scr
set MOVIES_SCRIP=%SCRDIR%\aux_movies.py
set FIGURES_SCRIPT=%SCRDIR%\aux_figures.py
set PATHWAYS_SCRIPT=%SCRDIR%\aux_pathways.py
set EXCITON_SCRIPT=%SCRDIR%\aux_excitons.py

rem use NUMBER_OF_PROCESSES to set PARALLEL
if %NUMBER_OF_PROCESSES% gtr 1 (
   set PARALLEL=True
) else (
   set PARALLEL=False
)

rem test MPI presence
set MPI_REPORT=MPI presence not tested
if %PARALLEL%==True (
   %PYTHON% %SCRDIR%\probe_mpi.py
   if %ERRORLEVEL% geq 1 (
      set MPI_REPORT=MPI not found: mpi4py package or MPI implementation missing
      set PARALLEL=False
   ) else (
      set MPI_REPORT=MPI probed with success!
      set PARALLEL=True
   )
)

rem set SAVE_OUTPUT depending on LOG_FILE

if [%LOG_FILE%] == [] (
   set SAVE_OUTPUT=False
) else (
   set SAVE_OUTPUT=True
)

if %PARALLEL% == True (
   set PARALLELOPT= -p -n %NUMBER_OF_PROCESSES% -d mpiexec
) else (
   set PARALLELOPT=
)

if %BACKGROUND% == True (
   set STARTER=start /b
) else (
   set STARTER=
)

if %SAVE_OUTPUT% == True (
   set PIPE= ^> %LOG_FILE%
) else (
   if %BACKGROUND% == True (
      set PIPE= ^> output.log
   ) else (
      set PIPE=
   )
)

rem ######################################################################
rem
rem     Tasks
rem
rem ######################################################################

rem     Default task
if "%1" == "" (
   set task=help
) else (
   set task=%1
)

rem     Running the simulation
if %task% == run (

   if %PARALLEL% == True (
      echo Running parallel simulation
   ) else (
      echo Running serial simulation
   )
   if %SAVE_OUTPUT% == True (
      echo Logging output into file: %LOG_FILE%
   )
   echo Issuing command: "%STARTER% qrhei run %PARALLELOPT% script_Policht2021.yaml %PIPE%"
   if %BACKGROUND% == True (
      echo Simulation runs in backgroud. Use taskmanager to watch its progress.
   ) else (
      echo Wait untill the simulation is finished ^(use Ctrl-c to stop execution^) ...
   )
   %STARTER% qrhei run %PARALLELOPT% %SCRPTNAME%.yaml %PIPE%


) else if %task% == figures (

   %PYTHON% %FIGURES_SCRIPT% %2

) else if %task% == movies (

   %PYTHON% %MOVIES_SCRIP% %2 %NUMBER_OF_PROCESSES%

) else if %task% == pathways (

   %PYTHON% %PATHWAYS_SCRIPT% %2 %3

) else if %task% == excitons (

   %PYTHON% %EXCITON_SCRIPT% %2

) else if %task% == copy (

   mkdir ..\%2
   copy  script_Followup2021.* ..\%2
   copy  make.bat ..\%2   
   copy  Makefile ..\%2
   copy  scr ..\%2
    


rem     Cleaning files
) else if %task% == clean (

   echo Removing files
   for /f %%i in ('dir /a:d /b sim_*') do rmdir /s /q %%i
   if exist *~ del *~
   if exist log del log
   if exist output.log del output.log
   if exist *.bak del *.bak

rem     Cleaning media files
) else if %task% == del (

   echo Removing files
   if exist *.png del *.png
   if exist *.mov del *.mov
   if exist *.mp4 del *.mp4
   if exist *.dat del *.dat
   if exist *.qrp del *.qrp
   if exist test.log del test.log

rem     Removing all files
) else if %task% == purge (
   make clean
   make del

rem     Backing up imput file
) else if %task% == back (
   if exist script_Followup2021.yaml move script_Followup2021.yaml script_Followup2021.yaml.bak

rem     Setting test in simulation mode single
) else if %task% == set_test_single (
   make back
   copy templates\script_Followup2021_test_single.yaml .\script_Followup2021.yaml

rem     Setting test in simulation mode disorder
) else if %task% == set_test_disorder (
   make back
   copy templates\data_test_disorder\random_state.qrp .\
   copy templates\script_Followup2021_test_disorder.yaml .\script_Followup2021.yaml

rem     Setting test in simulation mode scan
) else if %task% == set_test_scan (
   make back
   copy templates\script_Followup2021_test_scan.yaml .\script_Followup2021.yaml

rem     Validating results against saved data
) else if %task% == validate (
   %PYTHON% %VALIDATION_SCRIPT% %2

rem     Test scan
) else if %task% == test_scan (
   make set_test_scan
   make run
   make validate
   make clean
   echo test_scan ended with success >> test.log

rem     Test single
) else if %task% == test_single (
   make set_test_single
   make run
   make validate
   make clean
   echo test_single ended with success >> test.log

rem     Test disorder
) else if %task% == test_disorder (
   make set_test_disorder
   make run
   make validate
   make clean
   echo test_disorder ended with success >> test.log

rem     Complete set of tests
) else if %task% == test (
   make purge
   make test_single
   make test_disorder
   make test_scan

rem     Setting example in simulation mode single
) else if %task% == set_example_single_dimer (
   make back
   copy templates\script_Followup2021_example_single_dimer.yaml .\script_Followup2021.yaml

rem     Setting example in simulation mode disorder
) else if %task% == set_example_single (
   make back
   copy templates\script_Followup2021_example_disorder.yaml .\script_Followup2021.yaml

rem     Setting example in simulation mode scan
) else if %task% == set_example_scan (
   make back
   copy templates\script_Followup2021_example_scan.yaml .\script_Followup2021.yaml

rem    Help message
) else if %task% == help (
   echo.
   echo Simulation Makefile
   echo ===================
   echo.
   echo To configure the session, edit the switches at the start
   echo of this Makefile. You can choose between serial and parallel
   echo simulation and adjust logging and whether to run in the
   echo background or not
   echo.
   echo Available tasks:
   echo ----------------
   echo.
   echo ^> make help
   echo.
   echo     Prints this message
   echo.
   echo ^> make run
   echo.
   echo     Runs the simulations. On Linux and macOS
   echo     the job submission can be configured
   echo.
   echo ^> make figures results_directory
   echo.
   echo     Produces 2D omega_2 map figures from the simulation
   echo     which ran in the simulation mode "single" or "disorder"
   echo     ^(see configureation yaml file^). results_directory is
   echo     the directory containing results of Quantarhei simulation.
   echo.
   echo ^> make movies results_directory
   echo.
   echo     Produces movies for the energy gap scan from the simulation
   echo     which ran in the simulation mode "scan"
   echo     ^(see configureation yaml file^). results_directory is
   echo     the directory containing results of Quantarhei simulation.
   echo.
   echo ^> make pathways pws_file n
   echo.
   echo     Prints out Feynman diagrams of rephasing pathways contributing to
   echo     a spectral region specified within the script file
   echo     'scr/aux_pathways.py'. Integer number n specifies the maximum number 
   echo     of printed pathways.  
   echo.
   echo ^> make clean
   echo.
   echo     Deletes the output of the simulations
   echo.
   echo ^> make del
   echo.
   echo     Deletes media files created by scripts
   echo.
   echo.
   echo Current settings:
   echo -----------------
   echo.
   echo NUMBER_OF_PROCESSES=%NUMBER_OF_PROCESSES%
   echo %MPI_REPORT%
   echo BACKGROUND=%BACKGROUND%
   echo LOG_FILE=%LOG_FILE%
   echo PYTHON=%PYTHON%
   echo.
   echo Will run in parallel: ^(PARALLEL=^) %PARALLEL%
   echo Output will be saved: ^(SAVE_OUTPUT=^) %SAVE_OUTPUT%
   echo.
)
rem ##################################################################################
rem   EOF

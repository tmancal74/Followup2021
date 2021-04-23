# -*- coding: utf-8 -*-
"""
    Intro
    -----

    Simulation script as a supporting information for the manuscript:

    Veronica R. Policht, Andrew Niedringhaus, Cameron Spitzfaden,
    Philip D. Laible, David F. Bocian, Christine Kirmaier, Dewey Holten,
    Tomas Mancal and Jennifer P. Ogilvie,
    Hidden Vibronic and Excitonic Structure and Vibronic Coherence Transfer
    in the Bacterial Reaction Center, submitted 2021

    This script performs a set of simulations reproducing the numerical results
    presented in the above manuscript. The full set of files needed to run the
    simulations consists of:

    script_Policht2021.py ..... Python script
    script_Policht2021.yaml ... input file

    Some additional Python files are required to transform the raw output into
    figures, movies and text data files.

    The script is based on Quantarhei package
    (see http://github.com/tmancal74/quantarhei) which requires Python 3
    installation and the standard Numpy and Scipy packages to support its
    numerics. Several smaller packages are needed for some specific tasks
    related to data analysis and display.


    Installation instructions
    -------------------------

    1. The script requires Python 3. We recommend to install Anaconda individual
    distribution of Python. This distribution can be obtained for Linux, Mac and
    Windows systems, and it contains most of the required libraries right after
    its installation. It enables quick and easy installation of Quantarhei
    package itself.

    See https://www.anaconda.com/products/individual for download and
    installation instructions.

    2. Once you set up your Anaconda Python, you can install Quantarhei package
    by typing the following into your terminal (the character ">" represents
    your terminal's prompt)

    > conda install -c tmancal74 -c conda-forge quantarhei

    The latest available version of Quantarhei, together with its dependencies
    will be installed.

    (If you have a different distribution of Python, you can alternatively use
    the "pip" tool to install Quantarhei as:

    > pip install quantarhei

    )

    3. You might need to install additional packages for some specific tasks.
    For instance, if you want to use multiple processors to speed up some of
    the calculations, you might need to install the mpi4py package by:

    > conda install mpi4py

    In addition you need to have message passing interface (MPI) implementation
    installed.On Windows, we recommand to install Microsoft's implementation of
    the MPI, which can be found at
    https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi

    To produce movies and figures out of the data, your system needs to have
    Latex installed.
    On Mac, look for MaxTex (https://www.tug.org/mactex/mactex-download.html),
    Windows systems work well with MiKTex (https://miktex.org/download).

    To produce movies, ffmpeg program is required. Although there is a package
    called ffmpeg in conda, on some platforms (on Mac and Windows) it seems to
    be easier to install ffmpeg independentently of Python. The program has to
    be in your system's path.

    4. Once Quantarhei package is installed, test its presence on your computer.
    If you type:

    > qrhei --version

    into your terminal, you should get a message similar to

    Quantarhei package version:  0.0.64

    This means that the Quantarhei package of version 0.0.64 is present. Version
    number 0.0.64 is the minimum requirement to run the scripts. Newer versions
    should give you the same results. Should the syntax of quantarhei commands
    change during future development, this script would likely fail. Try to get
    an updated version by typing the following on the command line:

    > qrhei fetch -e script_Policht2021

    This should give you the version of the script compatible with the latest
    Quantarhei development. If even this fails, report a bug at

    https://github.com/tmancal74/quantarhei/issues

    or write to the author at tmancal74@gmail.com.


    Running the simulations
    -----------------------

    The simulation can be fully controlled by editing the input file

    script_Policht2021.yaml

    which contains both the configuration parameters of the simulation, and
    the numerical parameters of the model.

    In general, the script can be run by typing the following command into your
    terminal:

    > qrhei run script_Policht2021.yaml

    PLEASE, NOTE THAT THE INPUT FILE IS USED AS A PARAMETER OF THE qrhei SCRIPT

    If you want to run the simulations in parallel (this makes sense e.g.
    when averaging over energetic disorder), use the following command:

    > qrhei run -p -n 8 script_Policht2021.yaml

    This will run 8 copies of the script all of which will communicate with each
    other and share data to create a single output. The option -p is required
    because it tells qrhei driver to run in parallel. The option -n specifies
    the number of processes that will simultaneously start on your computer.

    To simplify running the code, we have provided Makefile to be run by the
    "make" utility (on Linux and Mac), and a batch file called make.bat to be
    run on Windows. The Windows version uses Makefile as a configuration file,
    so that configuration of the simulation is the same on all platforms. See
    Makefile for more details. One you decide whether you want to run in
    parallel or not and similar other decisiosn you make in Makefile, you just
    type:

    > make run

    to run the simulation. On Windows you can click runme.bat file in the
    simulation directory in stead.

    Data analysis can be also run by make.


    Data analysis and production of figures
    ---------------------------------------

    The simulation script produces a single output directory which contains
    raw data of the simulation. Depending on the use case, additional scripts
    are used to produce figures and movies representing the results. We provide
    two script to be found in the "aux" subdirectory.

    Typing on Linux/Max

    > make figures DIR=results_directory

    where results_directory is the name of the directory with Quantarhei
    output, creates figures (\omega_2 maps) corresponding to the simulation
    results.

    On Windows, equivalent command is simpler and reads

    > make figures results_directory

    Some simulations (the energy gap scans) can be presented by movies. The
    movies can be created from the Quantarhei output by typing:

    > make movies DIR=results_directory

    on Linux/Mac, and

    > make movies results_directory

    on Windows.


"""

################################################################################
################################################################################
#
#  PYTHON SCRIPT STARTS HERE
#
################################################################################
################################################################################

#
# IMPORTS
#

# standard Python libraries
import time
import datetime
import os
import shutil
import gc
import platform

# Numpy library
import numpy
import matplotlib.pyplot as plt

# Quantarhei imports
import quantarhei as qr
from quantarhei.utils.vectors import X
import quantarhei.functions as func
from quantarhei.core.units import kB_int
from quantarhei import printlog as print

print("\n*****   RC Simulation Script   *****")
print("\nUsing Quantarhei version", qr.Manager().version)

#
# SCRIPT INPUT FILE NAME
#
input_file = "input.yaml"
INP = qr.Input(input_file, show_input=False)

################################################################################
################################################################################
#
# SOME DEFINITIONS AND HELPER FUNCTIONS
#
################################################################################
################################################################################

def init_containers():
    """Initialization of the spectra containers

    """
    cont_p_re = qr.TwoDSpectrumContainer()
    cont_p_re.use_indexing_type("integer")
    cont_p_nr = qr.TwoDSpectrumContainer()
    cont_p_nr.use_indexing_type("integer")
    cont_m_re = qr.TwoDSpectrumContainer()
    cont_m_re.use_indexing_type("integer")
    cont_m_nr = qr.TwoDSpectrumContainer()
    cont_m_nr.use_indexing_type("integer")

    return (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr)


def save_containers(cont, dname, node=0):
    """Saves the content of containers

    """
    (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr) = cont
    name1 = "cont_p_re_"+str(node)
    name2 = "cont_p_nr_"+str(node)
    name3 = "cont_m_re_"+str(node)
    name4 = "cont_m_nr_"+str(node)
    fname = os.path.join(dname, name1)
    cont_p_re.savedir(fname)
    fname = os.path.join(dname, name2)
    cont_p_nr.savedir(fname)
    fname = os.path.join(dname, name3)
    cont_m_re.savedir(fname)
    fname = os.path.join(dname, name4)
    cont_m_nr.savedir(fname)

def save_averages(cont, dname):
    """Saves disorder averaged spectra

    """
    (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr) = cont
    name1 = "ave_p_re.qrp"
    name2 = "ave_p_nr.qrp"
    name3 = "ave_m_re.qrp"
    name4 = "ave_m_nr.qrp"
    fname = os.path.join(dname, name1)
    cont_p_re.save(fname)
    fname = os.path.join(dname, name2)
    cont_p_nr.save(fname)
    fname = os.path.join(dname, name3)
    cont_m_re.save(fname)
    fname = os.path.join(dname, name4)
    cont_m_nr.save(fname)


def unite_containers(node=0):
    """Collects container parts from the directory to make a single container

    """
    (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr) = init_containers()
    name1 = "cont_p_re_"+str(node)
    name2 = "cont_p_nr_"+str(node)
    name3 = "cont_m_re_"+str(node)
    name4 = "cont_m_nr_"+str(node)
    fname1 = os.path.join(dname, name1)
    cont_p_re = cont_p_re.unitedir(fname1)
    cont_p_re.save(fname1+".qrp")
    fname2 = os.path.join(dname, name2)
    cont_p_nr = cont_p_nr.unitedir(fname2)
    cont_p_nr.save(fname2+".qrp")
    fname3 = os.path.join(dname, name3)
    cont_m_re = cont_m_re.unitedir(fname3)
    cont_m_re.save(fname3+".qrp")
    fname4 = os.path.join(dname, name4)
    cont_m_nr = cont_m_nr.unitedir(fname4)
    cont_m_nr.save(fname4+".qrp")


def save_spectra(cont, ext="dat", dir=".", stype=qr.signal_TOTL):
    # saving total spectra
    drnm = os.path.join(dir, "spectra")
    try:
        os.makedirs(drnm)
    except FileExistsError:
        # directory already exists
        pass
    scont = cont.get_TwoDSpectrumContainer(stype=stype)
    tags = scont.tags
    for tg in tags:
        #sp.plot(show=True)
        sp = scont.get_spectrum(tag=tg)
        flnm = os.path.join(drnm, "sp_"+str(tg)+"."+ext)
        fgrn = os.path.join(drnm, "sp_"+str(tg)+".png")
        sp.plot(show=False)
        sp.savefig(fgrn)
        plt.close()
        print("Saving "+flnm)
        sp.save_data(flnm)
        if ext == "dat":
            _data = numpy.loadtxt(flnm, dtype=complex)
            print("max=", numpy.max(_data))

    try:
        save_cont = INP.total_spectrum["native_container"]
    except:
        save_cont = False

    if save_cont:
        cont.save(os.path.join(dir, "spectra_container.qrp"))




#
################################################################################
################################################################################
#
#  Main simulation routine
#
################################################################################
################################################################################
#
def run(omega, HR, dE, JJ, rate, E0, vib_loc="up", use_vib=True,
        detailed_balance=False, temperature=77.0, stype=qr.signal_REPH,
        save_eUt=False, t2_save_pathways=[], dname=None, trimer=None,
        disE=None):
    """Runs a complete set of simulations for a single set of parameters


    If disE is not None it tries to run averaging over Gaussian energetic
    disorder.

    """
    if dname is None:
        dname = "sim_"+vib_loc

    use_trimer =  trimer["useit"]
    rate_sp = trimer["rate"]

    #
    #  PARAMETERS FROM INPUT FILE
    #
    dip1 = INP.dip1 # [1.5, 0.0, 0.0]
    dip2 = INP.dip2 # [-1.0, -1.0, 0.0]
    width = INP.feature_width # 100.0
    width2 = INP.feature_width2

    normalize_maps_to_maximu = False
    trim_maps = False

    units = "1/cm"
    with qr.energy_units(units):

        data_descr = "_dO="+str(dE)+"_omega="+str(omega)+ \
                     "_HR="+str(HR)+"_J="+str(JJ)

        if use_vib:
            sys_char = "_vib"
        else:
            sys_char = "_ele"
        data_ext = sys_char+".png"
        obj_ext = sys_char+".qrp"

    # parameters of the SP
    if use_trimer:
        E2 = trimer["E2"]
        epsa = (E0+E2)/2.0
        DE = trimer["DE"]
        J2 = 0.5*numpy.sqrt(((E0-E2)**2)-(DE**2))
        ESP2 = epsa + DE/2.0
        ESP1 = epsa - DE/2.0

    #
    #   Model system is a dimer of molecules
    #
    with qr.energy_units("1/cm"):

        if not use_trimer:

            if disE is not None:
                mol1 = qr.Molecule([0.0, E0+disE[0]])
                mol2 = qr.Molecule([0.0, E0+dE+disE[1]])
                print("Monomer 1 energy:", E0+disE[0], "1/cm")
                print("Monomer 2 energy:", E0+dE+disE[1], "1/cm")
            else:
                mol1 = qr.Molecule([0.0, E0])
                mol2 = qr.Molecule([0.0, E0+dE])
                print("Monomer 1 energy:", E0, "1/cm")
                print("Monomer 2 energy:", E0+dE, "1/cm")
        else:

            if disE is not None:
                mol1 = qr.Molecule([0.0, ESP2+disE[0]])
                mol2 = qr.Molecule([0.0, E0+dE+disE[1]])
                print("Monomer 1 (SP_high) energy:", ESP2+disE[0], "1/cm")
                print("Monomer 2 (B) energy:", E0+dE+disE[1], "1/cm")
            else:
                mol1 = qr.Molecule([0.0, ESP2])
                mol2 = qr.Molecule([0.0, E0+dE])
                print("Monomer 1 (SP_high) energy:", ESP2, "1/cm")
                print("Monomer 2 (B) energy:", E0+dE, "1/cm")
            if disE is not None:
                mol3 = qr.Molecule([0.0, ESP1+disE[2]])
                print("Monomer 3 (SP_low) energy:", ESP1+disE[2], "1/cm")
            else:
                mol3 = qr.Molecule([0.0, ESP1])
                print("Monomer 3 (SP_low) energy:", ESP1, "1/cm")
            mol3.set_transition_width((0,1), width2)
            mol3.set_dipole(0,1, trimer["dipsp"])

        mol1.set_transition_width((0,1), width2)
        mol1.set_dipole(0,1, dip1)

        mol2.set_transition_width((0,1), width)
        mol2.set_dipole(0,1, dip2)

    if use_trimer:
        agg = qr.Aggregate([mol1, mol2, mol3])
    else:
        agg = qr.Aggregate([mol1, mol2])

    if use_trimer:

        with qr.energy_units("1/cm"):
            agg.set_resonance_coupling(0,1,JJ)
            print("B - SP_high coupling:", JJ, "1/cm")
            agg.set_resonance_coupling(0,2,J2)
            print("SP coupling:", J2, "1/cm")

    else:

        with qr.energy_units("1/cm"):
            agg.set_resonance_coupling(0,1,JJ)

    #
    # Electronic only aggregate
    #
    agg_el = agg.deepcopy()

    #
    # if nuclear vibrations are to be added, do it here
    #
    if use_vib:

        with qr.energy_units("1/cm"):
            mod1 = qr.Mode(omega)
            mod2 = qr.Mode(omega)

        if vib_loc == "down":
            set_vib = [True, False]
        elif vib_loc == "up":
            set_vib = [False, True]
        elif vib_loc == "both":
            set_vib = [True, True]
        else:
            raise Exception("Unknown location of the vibrations")

        if set_vib[0]:
            print("Vibrations set for SP_high molecule")
            mol1.add_Mode(mod1)
            mod1.set_nmax(0, INP.vibmode["no_g_vib"])
            mod1.set_nmax(1, INP.vibmode["no_e_vib"])
            mod1.set_HR(1, HR)

        if set_vib[1]:
            print("Vibrations set for B molecule")
            mol2.add_Mode(mod2)
            mod2.set_nmax(0, INP.vibmode["no_g_vib"])
            mod2.set_nmax(1, INP.vibmode["no_e_vib"])
            mod2.set_HR(1, HR)

    #
    # Before we build the aggregate, we make its copy to have an unbuilt version
    #
    agg3 = agg.deepcopy()
    aggA = agg.deepcopy()

    #
    # here we build the complete aggregate and its electronic only version
    #
    agg.build(mult=1)
    agg_el.build(mult=1)
    aggA.build(mult=1)


    # total Hamiltonian
    HH = agg.get_Hamiltonian()
    # electronic Hamiltonian
    He = agg_el.get_Hamiltonian()

    #
    # Laboratory setup
    #
    lab = qr.LabSetup()
    lab.set_polarizations(pulse_polarizations=[X,X,X],
                          detection_polarization=X)

    #
    # Time axis for the calculation of excited state evolution
    #
    t2_N_steps = INP.t2_N_steps
    t2_time_step = INP.t2_time_step
    time2 = qr.TimeAxis(0.0, t2_N_steps, t2_time_step)

    #
    # Containers for 2D maps with positive and negative frequencies
    #
    cont_p = qr.TwoDResponseContainer(t2axis=time2)
    cont_m = qr.TwoDResponseContainer(t2axis=time2)
    cont_tot = qr.TwoDResponseContainer(t2axis=time2)

    #
    # spectra will be indexed by the times in the time axis `time2`
    #
    cont_p.use_indexing_type(time2)
    cont_m.use_indexing_type(time2)
    cont_tot.use_indexing_type(time2)

    #
    # We define two-time axes, which will be FFTed and will define
    # the omega_1 and omega_3 axes of the 2D spectrum
    #
    t1_N_steps = INP.t1_N_steps
    t1_time_step = INP.t1_time_step
    t3_N_steps = INP.t3_N_steps
    t3_time_step = INP.t3_time_step
    t1axis = qr.TimeAxis(0.0, t1_N_steps, t1_time_step)
    t3axis = qr.TimeAxis(0.0, t3_N_steps, t3_time_step)

    #
    # This calculator calculated 2D spectra from the effective width
    #
    msc = qr.MockTwoDResponseCalculator(t1axis, time2, t3axis)
    with qr.energy_units("1/cm"):
        msc.bootstrap(rwa=E0, shape="Gaussian")

    #
    # System-bath interaction including vibrational states
    #
    operators = []
    rates = []

    print("Relaxation rates: ", rate, rate_sp, "1/fs")

    with qr.eigenbasis_of(He):

        if use_trimer:
            if He.data[3,3] < He.data[2,2]:
                Exception("Electronic states not orderred!")
            operators.append(qr.qm.ProjectionOperator(2, 3, dim=He.dim))
            with qr.energy_units("1/cm"):
                print("2 <- 3 : energies = ", He.data[2,2], He.data[3,3], "1/cm")
            rates.append(rate)
            print("Transfer time 2 <- 3:", 1.0/rate, "fs")
            
            if He.data[2,2] < He.data[1,1]:
                Exception("Electronic states not orderred!")
            operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
            with qr.energy_units("1/cm"):
                print("1 <- 2 : energies = ", He.data[1,1], He.data[2,2], "1/cm")
            rates.append(rate_sp)
            print("Transfer time 1 <- 2", 1.0/rate_sp, "fs")

        else:

            if He.data[2,2] < He.data[1,1]:
                Exception("Electronic states not orderred!")
            operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
            with qr.energy_units("1/cm"):
                print("1 <- 2 : energies = ", He.data[1,1], He.data[2,2], "1/cm")
            rates.append(rate_sp)
            print("Transfer time 1 <- 2", 1.0/rate, "fs")
            

    # include detailed balace
    if detailed_balance:
        if use_trimer:
            with qr.eigenbasis_of(He):
                T = temperature #77.0
                Den = (He.data[3,3] - He.data[2,2])/(kB_int*T)
                operators.append(qr.qm.ProjectionOperator(3, 2, dim=He.dim))
                thermal_fac = numpy.exp(-Den)
                rates.append(rate*thermal_fac)
                Den = (He.data[2,2] - He.data[1,1])/(kB_int*T)
                operators.append(qr.qm.ProjectionOperator(2, 1, dim=He.dim))
                thermal_fac = numpy.exp(-Den)
                rates.append(rate_sp*thermal_fac)

        else:

            with qr.eigenbasis_of(He):
                T = temperature #77.0
                Den = (He.data[2,2] - He.data[1,1])/(kB_int*T)
                operators.append(qr.qm.ProjectionOperator(2, 1, dim=He.dim))
                thermal_fac = numpy.exp(-Den)
            rates.append(rate*thermal_fac)

    sbi = qr.qm.SystemBathInteraction(sys_operators=operators, rates=rates)
    sbi.set_system(agg)

    #
    # Lindblad form for relaxation
    #
    LF = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)

    #
    # Pure dephasing
    #
    p_deph = qr.qm.ElectronicPureDephasing(agg, dtype="Gaussian")

    # we simplify calculations by converting dephasing to
    # corresponding Lorentzian form
    p_deph.convert_to("Lorentzian")

    eUt = qr.qm.EvolutionSuperOperator(time2, HH, relt=LF, pdeph=p_deph,
                                       mode="all")
    eUt.set_dense_dt(INP.fine_splitting)

    print("---")

    #
    # We calculate evolution superoperator
    #
    eUt.calculate(show_progress=False)

    # save the evolution operator
    if save_eUt:
        eut_name = os.path.join(dname, "eUt"+
                                    "_omega2="+str(omega)+data_descr+obj_ext)
        eUt.save(eut_name)

    #
    # Prepare aggregate with all states (including 2-EX band)
    #
    agg3.build(mult=2)
    agg3.diagonalize()

    agg_el.diagonalize()
    print("")
    print("Square of the transition dipoles")
    print(agg_el.D2)
    print("")
    print("---")


    #
    # calculation of absorption spectrum
    #
    time1 = qr.TimeAxis(0.0, 1000, 5.0)
    absc = qr.MockAbsSpectrumCalculator(time1, system=aggA)
    with qr.energy_units("1/cm"):
        absc.bootstrap(rwa=E0)

    spctrm = absc.calculate()
    spctrm.normalize2()

    with qr.energy_units("1/cm"):
        spctrm.plot(show=False, axis=[10500.0, 13500.0, 0.0, 1.1])
        spctrm.savefig(os.path.join(dname, "abs.png"))
        spctrm.save_data(os.path.join(dname, "abs.dat"))

    pways = dict()

    olow_cm = omega-INP.omega_uncertainty/2.0
    ohigh_cm = omega+INP.omega_uncertainty/2.0
    olow = qr.convert(olow_cm, "1/cm", "int")
    ohigh = qr.convert(ohigh_cm, "1/cm", "int")

    for t2 in time2.data:

        # this could save some memory of pathways become too big
        pways = dict()

        print("T2 =", t2, "fs (of T2_max =", time2.max, "fs)")

        twod = msc.calculate_one_system(t2, agg3, eUt, lab, pways=pways,
                                        dtol=1.0e-12,
                                        selection=[["omega2",[olow, ohigh]]])
        pws = pways[str(t2)]
        npa = len(pws)
        has_R = False
        has_NR = False
        for pw in pws:
            if pw.pathway_type == "NR":
                has_NR = True
            elif pw.pathway_type == "R":
                has_R = True

        if t2 in t2_save_pathways:
            pws_name = os.path.join(dname, "pws_t2="+str(t2)+
                                    "_omega2="+str(omega)+data_descr+obj_ext)
            qr.save_parcel(pways[str(t2)], pws_name)

        cont_p.set_spectrum(twod)

        twod = msc.calculate_one_system(t2, agg3, eUt, lab, pways=pways,
                                        dtol=1.0e-12,
                                        selection=[["omega2",[-ohigh, -olow]]])

        pws = pways[str(t2)]
        npa = len(pws)
        #print(" m:", npa)
        has_R = False
        has_NR = False
        for pw in pws:
            if pw.pathway_type == "NR":
                has_NR = True
            elif pw.pathway_type == "R":
                has_R = True

        #print(" R:", has_R, ", NR:", has_NR)

        if t2 in t2_save_pathways:
            pws_name = os.path.join(dname, "pws_t2="+str(t2)+
                                    "_omega2="+str(-omega)+data_descr+obj_ext)
            qr.save_parcel(pways[str(t2)], pws_name)

        cont_m.set_spectrum(twod)

        # calculation without pre-selecting pathways
        twod = msc.calculate_one_system(t2, agg3, eUt, lab, pways=pways,
                                        dtol=1.0e-12)

        cont_tot.set_spectrum(twod)

    # saving total spectrum to a directory for further analysis
    try:
        saveit = INP.total_spectrum["save_it"]
    except:
        saveit = False

    if saveit:
        try:
            dform = INP.total_spectrum["data_format"]
        except:
            dform = "dat"
        try:
            stp = INP.total_spectrum["spectra_type"]
            if stp == "TOTL":
                stype = qr.signal_TOTL
            elif stp == "REPH":
                stype = qr.signal_REPH
            elif stp == "NONR":
                stype = qr.signal_NONR
            else:
                raise Exception("Wrong signal type")
        except:
            stype = qr.signal_REPH

        save_spectra(cont_tot, ext=dform, dir=dname, stype=stype)


    #
    # Save aggregate when a single calculation is done
    #
    if save_eUt:
        fname = os.path.join(dname, "aggregate.qrp")
        agg3.save(fname)

    #
    # Window function for subsequenty FFT
    #
    window = func.Tukey(time2, r=INP.tukey_window_r, sym=False)

    #
    # FFT with the window function
    #
    # Specify REPH, NONR or `total` to get different types of spectra
    #
    print("Calculating FFT of the 2D maps")
    #fcont = cont.fft(window=window, dtype=stype) #, dpart="real", offset=0.0)

    fcont_p_re = cont_p.fft(window=window, dtype=qr.signal_REPH)
    fcont_p_nr = cont_p.fft(window=window, dtype=qr.signal_NONR)
    fcont_p_to = cont_p.fft(window=window, dtype=qr.signal_TOTL)

    if normalize_maps_to_maximu:
        fcont_p_re.normalize2(dpart=qr.part_ABS)
        fcont_p_nr.normalize2(dpart=qr.part_ABS)
        fcont_p_to.normalize2(dpart=qr.part_ABS)

    fcont_m_re = cont_m.fft(window=window, dtype=qr.signal_REPH)
    fcont_m_nr = cont_m.fft(window=window, dtype=qr.signal_NONR)
    fcont_m_to = cont_m.fft(window=window, dtype=qr.signal_TOTL)

    if normalize_maps_to_maximu:
        fcont_m_re.normalize2(dpart=qr.part_ABS)
        fcont_m_nr.normalize2(dpart=qr.part_ABS)
        fcont_m_to.normalize2(dpart=qr.part_ABS)

    if trim_maps:
        twin = INP.trim_maps_to
        with qr.energy_units("1/cm"):
            fcont_p_re.trimall_to(window=twin)
            fcont_p_nr.trimall_to(window=twin)
            fcont_p_to.trimall_to(window=twin)

    show_omega = omega

    with qr.frequency_units("1/cm"):
        sp1_p_re, show_Npoint1 = fcont_p_re.get_nearest(show_omega)
        sp2_p_re, show_Npoint2 = fcont_p_re.get_nearest(-show_omega)
        sp1_p_nr, show_Npoint1 = fcont_p_nr.get_nearest(show_omega)
        sp2_p_nr, show_Npoint2 = fcont_p_nr.get_nearest(-show_omega)
        sp1_p_to, show_Npoint1 = fcont_p_to.get_nearest(show_omega)
        sp2_p_to, show_Npoint2 = fcont_p_to.get_nearest(-show_omega)
        sp1_m_re, show_Npoint1 = fcont_m_re.get_nearest(show_omega)
        sp2_m_re, show_Npoint2 = fcont_m_re.get_nearest(-show_omega)
        sp1_m_nr, show_Npoint1 = fcont_m_nr.get_nearest(show_omega)
        sp2_m_nr, show_Npoint2 = fcont_m_nr.get_nearest(-show_omega)
        sp1_m_to, show_Npoint1 = fcont_m_to.get_nearest(show_omega)
        sp2_m_to, show_Npoint2 = fcont_m_to.get_nearest(-show_omega)

    sstm = platform.system()
    #print(sstm)
    if sstm != "Windows":
        import resource
        memo = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)
        print("Memory usage: ", memo, "in MB" )

    return (sp1_p_re, sp1_p_nr, sp2_m_re, sp2_m_nr)

#
#  END OF THE MAIN SIMULATION ROUTINE
#

###############################################################################
###############################################################################
#
#   START OF THE SCRIPT
#
###############################################################################
###############################################################################

#
# First collecting parameters of the simulation
#

# vibrational mode properties
parms1 = [INP.vibmode]
vib_loc = INP.location_of_vibrations

# this is a fix to have rate defined separately from other parameters
parms1[0]["rate"] = INP.rate

#
#   MODELS
#   This is meant to enable calculations of models with vibrations on different
#   molecules within one simulation. In the present script, we calculate only
#   a single model
#
models = [dict(vib_loc=vib_loc)]

#
# t2s at which pathways will be saved
#
t2_save_pathways = INP.t2_save_pathways #[50.0, 100.0, 200.0, 300.0]

#
# Here we construct a path through parameters space
#
center = INP.dE01 # This is the center of the disorder distribution
                  # of the energy gap between the reference monomer and the
                  # secondary monomer; or B and upper SP exciton in the case
                  # of the trimer model
step = INP.step #2.0

max_available_fwhm = INP.max_available_fwhm #100.0
how_many_fwhm = INP.how_many_fwhm #2

#step = 2
Ns_d = int(2.0*how_many_fwhm*max_available_fwhm/step) # 50
Ns_u = int(2.0*how_many_fwhm*max_available_fwhm/step) # 50

vax = qr.ValueAxis(center-Ns_d*step, Ns_d+Ns_u+1, step)
trimer = INP.trimer
use_trimer =  trimer["useit"]

#
# Here we specify pairs of parameters (resonance coupling J and energy
# gap \Delta E between the monomers). One could specify an arbitrary
# "pathway" in the parameters space. Below we specify a line of
# increasing \Delta E with constant J.
#
ptns = []

single_run = INP.single_realization
disorder = INP.disorder
detailed_balance = INP.detailed_balance
temperature = INP.temperature

#
# Run with a single realization (sigle set of parameters)
#
if single_run:

    ptns.append((INP.resonance_coupling, center, INP.trimer))

#
# Run with disorder and explicite averaging (sigle set + variations by disorder)
#
elif disorder:

    ptns.append((INP.resonance_coupling, center, INP.trimer))

#
# Many runs with predefined sets of parameters (for later averaging)
#
else:

    if use_trimer:

        for val in vax.data:
            ptns.append((INP.resonance_coupling, val,
                         INP.trimer))

    else:
        for val in vax.data:
            ptns.append((INP.resonance_coupling, val,
                         INP.trimer))

E0 = INP.E0 # transition energy (in 1/cm) of the reference monomer

#
# Containers for resulting 2D maps
#
if not disorder:
    (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr) = init_containers()

#
################################################################################
################################################################################
#
#   LOOP OVER SIMULATIONS
#
################################################################################
################################################################################
#

#
# Formally starting the regions where parallel loops and functions can be used
#
qr.start_parallel_region()
config = qr.distributed_configuration()

parms = parms1
i_p_re = 0
n_save = 0
tags = []
save_it_at_the_end = False

#
# loop over different models (we do not use it here - the goes only once)
#

tA = time.time()
at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
print("\nStarting simulation set at", at)
ll = 1
for model in models:
    print("Model no.", ll, "of", len(models))
    vib_loc = model["vib_loc"]
    try:
        app = INP.append_to_dirname
    except:
        app = ""
    dname = "sim_"+vib_loc+app
    ts = time.time()
    if INP.append_time_stamp:
        at = '_{0:%Y-%m-%d_%H%M%S}'.format(datetime.datetime.now())
        dname = dname+at
    try:
        os.makedirs(dname)
    except FileExistsError:
        # directory already exists
        pass

    if INP.copy_input_file_to_results:
        if INP._from_file:
            shutil.copy2(input_file, dname)
        else:
            INP.dump_yaml(os.path.join(dname, "input_file.yml"))

    #
    # Main loop
    #

    kk = 1
    at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    print("\nStarting the simulation at:", at)
    Np = len(parms)
    for par in parms:

        print("Run no.", kk, "of", Np)
        omega = par["omega"]
        HR = par["HR"]
        use_vib = par["use_vib"]
        rate = par["rate"]
        kp = 1
        Nje = len(ptns)

        #
        # loop over disorder
        #
        if disorder:

            (JJ, dE, trimer) = ptns[0]

            qr.timeit(show_stamp=True)

            Nreal = INP.N_realizations
            data_initialized = False
            disM = numpy.zeros((3,Nreal))
            sigma = INP.disorder_fwhm/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))

            if INP.random_state["reset"]:
                try:
                    random_state = qr.load_parcel(INP.random_state["file"])
                    numpy.random.set_state(random_state)
                except:
                    raise Exception("Loading random state failed")
            if INP.random_state["save"]:
                random_state = numpy.random.get_state()
                qr.save_parcel(random_state, INP.random_state["file"])

            for ri in range(Nreal):
                disM[:,ri] = sigma*numpy.random.randn(3)

            #
            # PARALLEL (if ON) LOOP OVER DISORDER
            #
            for ds in qr.block_distributed_range(0,Nreal):
                # generating random numbers
                disE = numpy.zeros(3,dtype=qr.REAL)

                if Nreal > 1:
                    disE[0] = disM[0,ds]
                    disE[1] = disM[1,ds]
                    disE[2] = disM[2,ds]

                print("\nCalculating disordered spectra ... (",ds+1,"of",Nreal,
                      ") [run ",kk,"of",Np,"]")
                print("---")
                print("Temperature =", temperature,"K")
                #print("JJ =", JJ, "1/cm")
                print("dE =", dE, "1/cm")
                print("Disorder in energies: ", disE, "1/cm")

                t1 = time.time()

                # only save the propagator if we calculate a single realization
                if (Nreal == 1):
                    save_eUt = True
                else:
                    save_eUt = False

                (sp1_p_re, sp1_p_nr, sp2_m_re, sp2_m_nr) = \
                run(omega, HR, dE, JJ, rate, E0, vib_loc, use_vib,
                    save_eUt=save_eUt,t2_save_pathways=t2_save_pathways,
                    dname=dname, trimer=trimer, disE=disE,
                    detailed_balance=detailed_balance,temperature=temperature)

                t2 = time.time()
                gc.collect()
                print("... done in",t2-t1,"sec")

                if not data_initialized:
                    params = dict(J=JJ, dE=dE, E0=E0, omega=omega,
                                  delta=INP.disorder_fwhm)
                    sp1_p_re.log_params(params)
                    sp1_p_nr.log_params(params)
                    sp2_m_re.log_params(params)
                    sp2_m_nr.log_params(params)

                    if INP.restart_disorder:
                        fname = os.path.join(dname, "ave_p_re.qrp")
                        av1_p_re = qr.load_parcel(fname)
                        fname = os.path.join(dname, "ave_p_nr.qrp")
                        av1_p_nr = qr.load_parcel(fname)
                        fname = os.path.join(dname, "ave_m_re.qrp")
                        av2_m_re = qr.load_parcel(fname)
                        fname = os.path.join(dname, "ave_m_nr.qrp")
                        av2_m_nr = qr.load_parcel(fname)
                    else:
                        av1_p_re = sp1_p_re.deepcopy()
                        av1_p_nr = sp1_p_nr.deepcopy()
                        av2_m_re = sp2_m_re.deepcopy()
                        av2_m_nr = sp2_m_nr.deepcopy()

                        av1_p_re.data[:,:] = 0.0
                        av1_p_nr.data[:,:] = 0.0
                        av2_m_re.data[:,:] = 0.0
                        av2_m_nr.data[:,:] = 0.0

                    data_initialized = True

                av1_p_re.data += sp1_p_re.data
                av1_p_nr.data += sp1_p_nr.data
                av2_m_re.data += sp2_m_re.data
                av2_m_nr.data += sp2_m_nr.data

                tags.append(i_p_re)

                i_p_re +=1
                kp += 1

            av1_p_re.data = config.reduce(av1_p_re.data)/Nreal
            av1_p_nr.data = config.reduce(av1_p_nr.data)/Nreal
            av2_m_re.data = config.reduce(av2_m_re.data)/Nreal
            av2_m_nr.data = config.reduce(av2_m_nr.data)/Nreal

            qr.finished_in(show_stamp=True)

        #
        # Loop over parameter sets
        #
        else:

            #
            # PARALLEL (if ON) LOOP OVER PARAMETER RANGE
            #
            for (JJ, dE, trimer) in qr.block_distributed_list(ptns):

                print("\nCalculating spectra ... (",kp,"of",Nje,
                      ") [run ",kk,"of",Np,"]")
                print("---")
                print("Temperature =", temperature,"K")
                #print("JJ =", JJ, "1/cm")
                print("dE =", dE, "1/cm")

                t1 = time.time()
                if Nje == 1:
                    save_eUt = True
                else:
                    save_eUt = False

                (sp1_p_re, sp1_p_nr, sp2_m_re, sp2_m_nr) = \
                run(omega, HR, dE, JJ, rate, E0, vib_loc, use_vib,
                    save_eUt=save_eUt, t2_save_pathways=t2_save_pathways,
                    dname=dname, trimer=trimer,
                    detailed_balance=detailed_balance, temperature=temperature)

                t2 = time.time()
                gc.collect()
                print("... done in",t2-t1,"sec")

                params = dict(J=JJ, dE=dE, E0=E0, omega=omega)
                sp1_p_re.log_params(params)
                sp1_p_nr.log_params(params)
                sp2_m_re.log_params(params)
                sp2_m_nr.log_params(params)

                cont_p_re.set_spectrum(sp1_p_re, tag=i_p_re)
                cont_p_nr.set_spectrum(sp1_p_nr, tag=i_p_re)
                cont_m_re.set_spectrum(sp2_m_re, tag=i_p_re)
                cont_m_nr.set_spectrum(sp2_m_nr, tag=i_p_re)
                tags.append(i_p_re)

                n_save += 1
                if not save_it_at_the_end:

                    if numpy.mod(n_save,10) == 0:
                        # we save and release containers after some time
                        print("Saving intermediate results; cleaning memory")
                        cont = (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr)
                        save_containers(cont, dname, node=config.rank)
                        (cont_p_re, cont_p_nr,
                        cont_m_re, cont_m_nr) = init_containers()

                i_p_re +=1
                kp += 1

        kk += 1
    ll += 1

tB = time.time()
at = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
print("\n... finished simulation set at", at, "in", tB-tA,"sec")

#
################################################################################
#  Final clean-up
################################################################################
#

if disorder:
    if config.rank == 0:
        cont = (av1_p_re, av1_p_nr, av2_m_re, av2_m_nr)
        save_averages(cont, dname)

else:
    if save_it_at_the_end:
        rank = config.rank
        fname = os.path.join(dname, "cont_p_re_"+str(rank)+".qrp")
        cont_p_re.save(fname)
        fname = os.path.join(dname, "cont_p_nr_"+str(rank)+".qrp")
        cont_p_nr.save(fname)
        fname = os.path.join(dname, "cont_m_re_"+str(rank)+".qrp")
        cont_m_re.save(fname)
        fname = os.path.join(dname, "cont_m_nr_"+str(rank)+".qrp")
        cont_m_nr.save(fname)

    else:
        # saving the rest of containers
        cont = (cont_p_re, cont_p_nr, cont_m_re, cont_m_nr)
        save_containers(cont, dname, node=config.rank)

        # uniting the containers saved in pieces into one file each
        unite_containers(node=config.rank)

#
# Formal closing of the region that can be run in parallel
#
qr.close_parallel_region()

################################################################################
################################################################################
# EOF

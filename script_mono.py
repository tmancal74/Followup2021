import platform
import os

import numpy
import quantarhei as qr
from quantarhei.utils.vectors import X
import quantarhei.functions as func


E0 = 10000.0
width=100.0

fine_splitting = 10
tukey_window_r = 0.3
normalize_maps_to_maximum = True
trim_maps = False
omega = 500.0
show_omega = omega

freqs = [500.0, 200.0, 300.0, 500.0, 700.0, 900.0]
hrs   = [0.5, 0.1, 0.1, 0.03,     0.03,   0.03]
Nmods = len(freqs)


Nmods = 2
Nmax_g = 3
Nmax_e = 3

with qr.energy_units("1/cm"):

    # two-level molecule
    mol = qr.Molecule([0.0, E0])
    mol.set_dipole(0,1,[1.0, 0.0, 0.0])

    # modes of vibratinal motion
    for ii in range(Nmods):
        frq = freqs[ii]
        hr = hrs[ii]

        mod = qr.Mode(frequency=frq)
        mol.add_Mode(mod)
        mod.set_HR(N=1, hr=hr)
        mod.set_nmax(N=1, nmax=Nmax_e)
        mod.set_nmax(N=0, nmax=Nmax_g)


# single member aggregate
# (needs to be created to allow spectroscopic calculations)
agg = qr.Aggregate(molecules=[mol])

# setting system bath interaction to provide lineshape
with qr.energy_units("1/cm"):
    mol.set_transition_width((0,1), width)

# building aggregate
agg.build()
HH = agg.get_Hamiltonian()

# check the aggregate
print(agg)

#
# calculation of absorption spectrum
#
time1 = qr.TimeAxis(0.0, 1000, 5.0)
absc = qr.MockAbsSpectrumCalculator(time1, system=agg)
with qr.energy_units("1/cm"):
    absc.bootstrap(rwa=E0)

spctrm = absc.calculate()
spctrm.normalize2()

with qr.energy_units("1/cm"):
    spctrm.plot(show=False, axis=[9000.0, 13000.0, 0.0, 1.1])
    spctrm.savefig("abs.png")

#
# calculation of 2D spectrum
#
time2 = qr.TimeAxis(0.0,10, 10.0)
time3 = qr.TimeAxis(0.0, time1.length, time1.step)

#
# Laboratory setup
#
lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X],
                        detection_polarization=X)

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

#
# This calculator calculated 2D spectra from the effective width
#
msc = qr.MockTwoDResponseCalculator(time1, time2, time3)
with qr.energy_units("1/cm"):
    msc.bootstrap(rwa=E0, shape="Gaussian")

#
# Pure dephasing
#
p_deph = qr.qm.ElectronicPureDephasing(agg, dtype="Gaussian")

# we simplify calculations by converting dephasing to
# corresponding Lorentzian form
p_deph.convert_to("Lorentzian")

operators=[]
operators.append(qr.qm.ProjectionOperator(2, 1, dim=2))

rates=[1.0/1000000.0]

sbi = qr.qm.SystemBathInteraction(sys_operators=operators, rates=rates)
sbi.set_system(agg)

#
# Lindblad form for relaxation
#
LF = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)

eUt = qr.qm.EvolutionSuperOperator(time2, HH, relt=LF, pdeph=p_deph,
                                    mode="all")
eUt.set_dense_dt(fine_splitting)

#
# We calculate evolution superoperator
#
eUt.calculate(show_progress=False)

olow_cm = omega-10.0/2.0
ohigh_cm = omega+10.0/2.0
olow = qr.convert(olow_cm, "1/cm", "int")
ohigh = qr.convert(ohigh_cm, "1/cm", "int")

print("---")
print("Calculating 2D spectra")
for t2 in time2.data:

    # this could save some memory of pathways become too big
    pways = dict()

    print("T2 =", t2, "fs (of T2_max =", time2.max, "fs)")

    twod = msc.calculate_one_system(t2, agg, eUt, lab, pways=pways,
                                    dtol=1.0e-12, selection=[["omega2",[olow, ohigh]]])
    pws = pways[str(t2)]
    npa = len(pws)
    #print(" p:", npa)
    has_R = False
    has_NR = False
    for pw in pws:
        if pw.pathway_type == "NR":
            has_NR = True
        elif pw.pathway_type == "R":
            has_R = True

    cont_p.set_spectrum(twod)


    twod = msc.calculate_one_system(t2, agg, eUt, lab, pways=pways,
                                    dtol=1.0e-12, selection=[["omega2",[-ohigh, -olow]]])

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

    cont_m.set_spectrum(twod)

    # calculate 2D spectra without pre-selecting pathways
    twod = msc.calculate_one_system(t2, agg, eUt, lab, pways=pways,
                                    dtol=1.0e-12)
    cont_tot.set_spectrum(twod)


def save_spectra(cont, ext="dat"):
    # saving total spectra
    drnm = "spectra"
    try:
        os.makedirs(drnm)
    except FileExistsError:
        # directory already exists
        pass
    scont = cont.get_TwoDSpectrumContainer()
    tags = scont.tags
    for tg in tags:
        sp = scont.get_spectrum(tag=tg)
        flnm = os.path.join(drnm, "sp_mono_"+str(tg)+"."+ext)
        fgrn = os.path.join(drnm, "sp_mono_"+str(tg)+".png")
        sp.plot(show=False)
        sp.savefig(fgrn)
        print("Saving "+flnm)
        sp.save_data(flnm)
        if ext == "dat":
            _data = numpy.loadtxt(flnm, dtype=complex)
            print("max=", numpy.max(_data))

save_spectra(cont_tot,"dat")


#
# Window function for subsequenty FFT
#
window = func.Tukey(time2, r=tukey_window_r, sym=False)

#
# FFT with the window function
#
# Specify REPH, NONR or `total` to get different types of spectra
#
print("Calculating FFT of the 2D maps")
#fcont = cont.fft(window=window, dtype=stype) #, dpart="real", offset=0.0)

print("Positive frequency:")
fcont_p_re = cont_p.fft(window=window, dtype=qr.signal_REPH)
print("# 1/3")
fcont_p_nr = cont_p.fft(window=window, dtype=qr.signal_NONR)
print("# 2/3")
fcont_p_to = cont_p.fft(window=window, dtype=qr.signal_TOTL)
print("# 3/3")

if normalize_maps_to_maximum:
    fcont_p_re.normalize2(dpart=qr.part_ABS)
    fcont_p_nr.normalize2(dpart=qr.part_ABS)
    fcont_p_to.normalize2(dpart=qr.part_ABS)

print("Negative frequency:")
fcont_m_re = cont_m.fft(window=window, dtype=qr.signal_REPH)
print("# 1/3")
fcont_m_nr = cont_m.fft(window=window, dtype=qr.signal_NONR)
print("# 2/3")
fcont_m_to = cont_m.fft(window=window, dtype=qr.signal_TOTL)
print("# 3/3")

if normalize_maps_to_maximum:
    fcont_m_re.normalize2(dpart=qr.part_ABS)
    fcont_m_nr.normalize2(dpart=qr.part_ABS)
    fcont_m_to.normalize2(dpart=qr.part_ABS)

if trim_maps:
    twin = INP.trim_maps_to
    with qr.energy_units("1/cm"):
        fcont_p_re.trimall_to(window=twin)
        fcont_p_nr.trimall_to(window=twin)
        fcont_p_to.trimall_to(window=twin)

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

with qr.frequency_units("1/cm"):
    window=[8000.0, 12000.0, 8000.0, 12000.0]
    sp1_p_re.plot(show=False, window=window)
    sp1_p_re.savefig("twod_p_re_plus.png")
    sp2_p_re.plot(show=False, window=window)
    sp2_p_re.savefig("twod_p_re_minus.png")
    sp1_p_nr.plot(show=False, window=window)
    sp1_p_nr.savefig("twod_p_nr_plus.png")
    sp2_p_nr.plot(show=False, window=window)
    sp2_p_nr.savefig("twod_p_nr_minus.png")

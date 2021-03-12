import numpy
import quantarhei as qr


width=100.0

freqs = [100.0, 200.0, 300.0, 500.0, 700.0, 900.0]
hrs   = [0.1, 0.1, 0.1, 0.03,     0.03,   0.03]
Nmods = len(freqs)


#Nmods = 1
Nmax = 3

with qr.energy_units("1/cm"):

    # two-level molecule
    mol = qr.Molecule([0.0, 10000.0])
    mol.set_dipole(0,1,[1.0, 0.0, 0.0])

    # modes of vibratinal motion
    for ii in range(Nmods):
        frq = freqs[ii]
        hr = hrs[ii]

        mod = qr.Mode(frequency=frq)
        mol.add_Mode(mod)
        mod.set_HR(N=1, hr=hr)
        mod.set_nmax(N=1, nmax=Nmax)


# single member aggregate 
# (needs to be created to allow spectroscopic calculations)
agg = qr.Aggregate(molecules=[mol])

# setting system bath interaction to provide lineshape
with qr.energy_units("1/cm"):
    mol.set_transition_width((0,1), width)

# building aggregate
agg.build()

# check the aggregate
print(agg)

# calculation of absorption spectrum
time1 = qr.TimeAxis(0.0, 1000, 5.0)
absc = qr.MockAbsSpectrumCalculator(time1, system=agg)
with qr.energy_units("1/cm"):
    absc.bootstrap(rwa=10000.0)

spctrm = absc.calculate()

with qr.energy_units("1/cm"):
    spctrm.plot(show=True, axis=[9000.0, 13000.0, 0.0, 15.0])





import sys
import numpy
import quantarhei as qr


ptype = "REPH"
window_L = [11500, 12500, 12000, 13000]
window_U = [11500, 12500, 13000, 13500]

window = window_U


fname = sys.argv[1]
try:
    N = int(sys.argv[2])
except:
    N = 5
    print("Using default value of N =", N)

pws = qr.load_parcel(fname)

Nparcel = len(pws)
print("Number of pathways loaded:", Nparcel)

with qr.energy_units("1/cm"):

    pan = qr.LiouvillePathwayAnalyzer(pathways=pws)
    pan.select_frequency_window(window=window,
            replace=True)
    pws = pan.select_type(ptype=ptype, replace=False)
    Nselect = len(pws)
print("Number of pathways selected:", Nselect)

Nshow = numpy.min([N, Nselect])
print("Showing:", Nshow)


with qr.energy_units("1/cm"):

    pan = qr.LiouvillePathwayAnalyzer(pathways=pws)
    pan.select_frequency_window(window=window,
            replace=True)
    pws = pan.select_type(ptype=ptype, replace=False)

    for ii in range(Nshow):
        print(pws[ii])
        #print(numpy.abs(pws[ii].pref))


    print("Absolute values of prefactors:")
    print("")
    for ii in range(Nshow):
        #print(pws[ii])

        apref = numpy.abs(pws[ii].pref)
        apref0 = numpy.abs(pws[0].pref)
        print(apref,apref/apref0)

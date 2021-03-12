
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import scipy.io as io
import numpy
import matplotlib as mpl
import sys
import os

import quantarhei as qr

################################################################################
#
#  ADVANCED CONFIGURATION SECTION
#
################################################################################

# Number of contours
Ncont = 10

# 2D spectrum range in 1/cm
window=[10500,13550,10500,13550]

#cmpfile = None
cmpfile = os.path.join("scr","parula_colormap.dat")

# movie format (mp4 or mov)
movie_ext = "mov"

###############################################################################
#
#  END OF CONFIGURATION SECTION
#
###############################################################################

try:
    target_dir = sys.argv[1]
except:
    print("Simulation output directory not specified")
    qr.exit()

try:
    Nnods = int(sys.argv[2])
except:
    print("Number of processes not specified")
    qr.exit()

ext = {0:"p_re", 1:"p_nr", 2:"m_re", 3:"m_nr"}
fig = None

# use submitted colormap
if cmpfile is not None:

    # conversion of matlab colormap to text
    #clrmp = io.loadmat(cmpfile)["colmap"]
    #cmap = mpl.colors.ListedColormap(clrmp)
    #numpy.savetxt("colormap.dat", clrmp)

    clrmp = numpy.loadtxt(cmpfile)
    cmap = mpl.colors.ListedColormap(clrmp)

else:
    cmap = None

do_nodes = True
normalize = True

def label_func(sp):
    """Function returning label for each frame
    """
    return (r'$\Delta E ='+str(sp.params["dE"])+r'\;\rm{cm}^{-1}$', [0.1,0.1])


def show_states_func(sp):
    """Function returning states for each frame
    """
    E0 = sp.params["E0"]
    dE = sp.params["dE"]
    omega = sp.params["omega"]
    return [[E0,"--k"], [E0+dE, "--b"], [E0+dE-omega, "--r"]]


if do_nodes:
    nodes = [i for i in range(Nnods)] #0, 1, 2, 3, 4]
else:
    nodes = [0]

print("\n*** Making energy gap scan movie ***\n")

for ext_i in ext:

    cont = qr.TwoDSpectrumContainer()
    cont.use_indexing_type("integer")


    #
    # Unite nodes
    #
    ii = 0
    print("Creating spectral container ...")
    for node in nodes:

        if do_nodes:
            ndp = "_"+str(node)
        else:
            ndp = ""

        file_name = os.path.join(target_dir, "cont_"+ext[ext_i]+ndp+".qrp")
        conta = qr.load_parcel(file_name)

        for tag in conta.spectra:
            #print(node, tag, ii)
            sp = conta.get_spectrum(tag)
            #print("dE = ", sp.params["dE"])
            ntag = ii

            if normalize:
                sp.normalize2(dpart=qr.part_ABS)

            cont.set_spectrum(sp, tag=ntag)
            ii += 1

    mfilename = "movie_"+ext[ext_i]+"_cont="+str(Ncont)+"."+movie_ext
    print("Exporting movie: ", mfilename)
    #
    # make a united movie
    #
    with qr.energy_units("1/cm"):
        cont.make_movie(mfilename, cmap=cmap,
        spart=qr.part_ABS, Npos_contours=Ncont, label_func=label_func,
        show_states_func=show_states_func, dpi=50)
    print("...done\n")

print("... finished")

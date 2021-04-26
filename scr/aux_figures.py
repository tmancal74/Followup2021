###############################################################################
#
#  This script produces figures from the raw Quantarhei data output
#
###############################################################################
#
#  HOW RUN THIS SCRIPT
#  -------------------
#
#  Use the Makefile in the same directory as this script. The usage
#  differs slightly between Windows and Linux/Unix/Mac systems.
#  On Linux/Unix/Mac systems you type
#
#  > make figures DIR=results_dir
#
#  the same on Windows reads
#
#  > make figures results_dir
#
#  In both cases results_dir is the name of the output directory of Quantarhei
#  simulation. In both cases this script will be run by the Python interpreter
#  and figures will be produced.
#
#
#  ADVANCED CONFIGURATION
#  ----------------------
#
#  In the advanced configuration section below you cn configure the output
#  figure
#
##  You can also supply a color map as a *.mat file (produced by Matlab) as
#  a second argument as below:
#
#  On Linux/Unix/Mac
#
#  > make figures DIR=results_dir CMAP=colormap_file
#
#  and
#
#  > make figures results_dir colormap_file
#
#  on Windows.
#

import scipy.io as io
import matplotlib as mpl
import sys, os
import numpy

import quantarhei as qr

###############################################################################
#
#  ADVANCED CONFIGURATION SECTION
#
###############################################################################

# Number of contours
Ncont = 10

# 2D spectrum range in 1/cm
window=[10500+500,13550+500,10500+500,13550+500]

points = [[12500, 13000],[13000,13000]]
maxs = [["square",[12400,12600,12800,13200]], ["square",[12800,13500,12800,13200]]]

#
# set lines along the x axis which should be plotted
# position along y axis in 1/cm
#
wlines=[13000.0, 13570.0]

#cmpfile = None
cmpfile = os.path.join("scr", "parula_colormap.dat") #"parula_colormap.mat"

###############################################################################
#
#  END OF CONFIGURATION SECTION
#
###############################################################################

print()
print("===============================")
print("      Creating figures")
print("===============================")

try:
    target_dir = sys.argv[1]
except:
    print("Simulation output directory not specified")
    qr.exit()

ext = {0:"p_re", 1:"p_nr", 2:"m_re", 3:"m_nr"}
fig = None

combinations = [["p", "re"], ["p", "nr"], ["m", "re"], ["m", "nr"]]

nodes = [0] #, 1, 2]

mode = "single"
if mode == "single":
    do_nodes = True
    prefix = "cont_"
    postfix = "_single"
elif mode == "average":
    do_nodes = False
    prefix = "ave_"
    postfix = "_average"

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

try:

    print("\nLooking for single realization files ...")
    for ext_i in range(4):

        print("\nSignal component:", ext[ext_i])

        cont = qr.TwoDSpectrumContainer()
        cont.use_indexing_type("integer")

        ii = 0
        for node in nodes:

            if do_nodes:
                ndp = "_"+str(node)
            else:
                ndp = ""

            file_name = os.path.join(target_dir, prefix+ext[ext_i]+ndp+".qrp")
            print("Loading file:", file_name)
            conta = qr.load_parcel(file_name)

            for tag in conta.spectra:
                #print("(node, tag, new tag):", node, tag, ii)
                sp = conta.get_spectrum(tag)
                ntag = ii
                cont.set_spectrum(sp, tag=ntag)
                ii += 1

        #print("Summary ("+ext[ext_i]+"):")
        for tag in cont.spectra:

            sp = cont.get_spectrum(tag)
            #print(tag, sp.params["dE"])
            
            with qr.energy_units("1/cm"):
                if fig is None:
                    fig = sp.plot(spart=qr.part_ABS, Npos_contours=Ncont,
                                  window=window,
                                  cmap=cmap, vmin_ratio=0.0)
                else:
                    fig  = sp.plot(fig=fig)

            file_name = ("fig_"+ext[ext_i]+postfix+"_cont="
                         +str(Ncont)+".png")
            print("Saving file:", file_name)
            sp.savefig(file_name)

            sval = []
            for area in maxs:
                with qr.energy_units("1/cm"):
                    pt = sp.get_area_max(area)
                    sval.append([area[1][0], area[1][1], area[1][2], area[1][3], numpy.abs(numpy.real(pt))])

            asval = numpy.array(sval)
            numpy.savetxt("max_"+ext[ext_i]+postfix+"_cont="
                         +str(Ncont)+".dat",asval)

            fig = None

#            sp.normalize2(dpart=qr.part_ABS)
#
#            with qr.energy_units("1/cm"):
#                if fig is None:
#                    fig = sp.plot(spart=qr.part_ABS, Npos_contours=Ncont,
#                                  window=window,
#                                  cmap=cmap, vmin_ratio=0.0)
#                else:
#                    fig  = sp.plot(fig=fig)
#
#            file_name = ("fig_"+ext[ext_i]+postfix+"_cont="
#                         +str(Ncont)+"_norm.png")
#            print("Saving file:", file_name)
#            sp.savefig(file_name)


            # taking cuts
            with qr.energy_units("1/cm"):

                for wline in wlines:
                    print("Spectral cut at "+str(wline))
                    cutfce = sp.get_cut_along_x(wline)
                    x = cutfce.axis
                    y = numpy.abs(numpy.real(cutfce.data))
                    cutfce._make_me(x,y)
                    cutfce.plot(show=False)
                    cutfce.savefig("cut_"+ext[ext_i]+"_"+str(wline)+".png")
                    cutfce.save_data("cut_"+ext[ext_i]+"_"+str(wline)+".dat",
                                     with_axis=cutfce.axis)


    print("\n... single realization files processed")
except:
    raise Exception()
    print("\nNo single spectrum files found")


try:

    print("\nLooking for averaged disorder files ...")
    dname = target_dir

    for comb in combinations:
        fsign = comb[0]
        tsigl = comb[1]

        fname = "ave_"+fsign+"_"+tsigl+".qrp"
        print("\nSignal component:", fsign+"_"+tsigl)

        floc = os.path.join(dname,fname)
        print("Loading file:", floc)
        av = qr.load_parcel(floc)

        mx = numpy.max(numpy.abs(av.data))
        av.data = av.data/mx

        try:
            with qr.energy_units("1/cm"):
                av.plot(spart=qr.part_ABS, Npos_contours=Ncont, window=window,
                cmap=cmap, vmin_ratio=0.0)
                flname = ("fig_"+fsign+"_"+tsigl+
                          "_average_cont="+str(Ncont)+".png")
                print("Saving file: ",flname)
                av.savefig(flname)
        except:
            raise Exception()

    print("\n... averaged spectra processed")
except:
    print("\nNo averaged spectrum files found")


print("\n       ... finished")
print("===============================")

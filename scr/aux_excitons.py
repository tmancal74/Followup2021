import sys
import os

import quantarhei as qr

try:
    dir = sys.argv[1]
except:
    print("Cannot find command line argument: Specify simulation directory.")
    qr.stop()

agpath = os.path.join(dir,"aggregate.qrp")
agg = qr.load_parcel(agpath)

agg.exciton_report(start=2)

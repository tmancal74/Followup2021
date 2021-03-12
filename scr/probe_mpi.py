"""
    Short script to probe presence of mpi4py package and mpi implementation
    
    
"""
import sys
import subprocess

return_value = 0

def _try_cmd(cmd):
    ret = False
    try:
        p = subprocess.Popen(cmd,
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
            
        # read and print output
        for line in iter(p.stdout.readline, b''):
            ln = line.decode()
            
        retval = p.wait()
    except:
        pass
    if retval == 0:
        ret = True
        
    return ret


# Test for mpi4py package
try:
    import mpi4py
except:
    print("Loading mpi4py failed!")
    print("Try 'conda install mpi4py' or 'pip install mpi4py' to install the package.") 
    return_value = 1
    
# Test for MPI's mpirun or mpiexec
if _try_cmd("mpirun --version") or _try_cmd("mpiexec"):
    return_value = 0
else:
    print("mpirun/mpiexec missing")
    return_value = 1

sys.exit(return_value)


    
    

import numpy as np, h5py, sys, getopt

def read_and_convert_T(tname):
   tfile = h5py.File(tname,"r")
   Taa_r = np.asarray(tfile['Taa_r'][:])
   Taa_i = np.asarray(tfile['Taa_i'][:])
   Tab_r = np.asarray(tfile['Tab_r'][:])
   Tab_i = np.asarray(tfile['Tab_i'][:])
   Tba_r = np.asarray(tfile['Tba_r'][:])
   Tba_i = np.asarray(tfile['Tba_i'][:])
   Tbb_r = np.asarray(tfile['Tbb_r'][:])
   Tbb_i = np.asarray(tfile['Tbb_i'][:])

def args(argv):
   tname = "T.h5"
   try:
      opts, args = getopt.getopt(argv,"i:")
   except getopt.GetoptError:
      print('convert_T-matrix.py -i <T-file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-i'):
         tname = arg
   return tname

if __name__ == "__main__":
   tname = args(sys.argv[1:])
   read_and_convert_T(tname)
      

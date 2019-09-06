import numpy as np, h5py, sys, getopt

ref_ind = 1.686 + 0.0312j

eps_r = np.real(ref_ind**2)
eps_i = np.imag(ref_ind**2)

def read_and_convert_T(tname,wavelen,a):
   tfile = h5py.File(tname+'.h5',"r")
   Taa_r = np.asarray(tfile['Taa_r'][:])
   nm = Taa_r.shape[0]
   Taai_r = np.reshape(Taa_r,nm**2,order='F')
   Taai_i = np.reshape(np.asarray(tfile['Taa_i'][:]),nm**2,order='F')
   Tabi_r = np.reshape(np.asarray(tfile['Tab_r'][:]),nm**2,order='F')
   Tabi_i = np.reshape(np.asarray(tfile['Tab_i'][:]),nm**2,order='F')
   Tbai_r = np.reshape(np.asarray(tfile['Tba_r'][:]),nm**2,order='F')
   Tbai_i = np.reshape(np.asarray(tfile['Tba_i'][:]),nm**2,order='F')
   Tbbi_r = np.reshape(np.asarray(tfile['Tbb_r'][:]),nm**2,order='F')
   Tbbi_i = np.reshape(np.asarray(tfile['Tbb_i'][:]),nm**2,order='F')
   
   with h5py.File(tname+"-fixed.h5","w") as f:
      dset1 = f.create_dataset("Taai_r", Taai_r.shape, dtype='double' )
      dset1[...] = Taai_r
      dset2 = f.create_dataset("Taai_i", Taai_i.shape, dtype='double' )
      dset2[...] = Taai_i
      dset3 = f.create_dataset("Tabi_r", Tabi_r.shape, dtype='double' )
      dset3[...] = Tabi_r
      dset4 = f.create_dataset("Tabi_i", Tabi_i.shape, dtype='double' )
      dset4[...] = Tabi_i
      dset5 = f.create_dataset("Tbai_r", Tbai_r.shape, dtype='double' )
      dset5[...] = Tbai_r
      dset6 = f.create_dataset("Tbai_i", Tbai_i.shape, dtype='double' )
      dset6[...] = Tbai_i
      dset7 = f.create_dataset("Tbbi_r", Tbbi_r.shape, dtype='double' )
      dset7[...] = Tbbi_r
      dset8 = f.create_dataset("Tbbi_i", Tbbi_i.shape, dtype='double' )
      dset8[...] = Tbbi_i     
      dset9 = f.create_dataset("T-ref-a", [3], dtype='double')
      dset9[...] = [eps_r, eps_i, a]
      dset10 = f.create_dataset("T-wavlens", [1], dtype='double')
      dset10[...] = wavelen

def args(argv):
   tname = "T"
   wavelen = 600e-9
   a = 1e-7
   try:
      opts, args = getopt.getopt(argv,"i:a:w:")
   except getopt.GetoptError:
      print('convert_T-matrix.py -i <T-file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-i'):
         tname = arg
      if opt in ('-a'):
         a = float(arg)
      if opt in ('-w'):
         wavelen = float(arg)
   return tname, wavelen, a

if __name__ == "__main__":
   tname, wavelen, a = args(sys.argv[1:])
   read_and_convert_T(tname,wavelen,a)
   
 

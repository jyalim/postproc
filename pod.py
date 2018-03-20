from hellaPy import *
from numpy import *
from pylab import *
from cheb import *
import glob,sys,natsort
import multiprocessing as mp
import ctypes
mkl_rt = ctypes.CDLL('libmkl_rt.so')
mkl_get_max_threads = mkl_rt.mkl_get_max_threads

def mkl_set_num_threads(cores):
    mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(cores)))

NPROC=50
SKIP=5
SKIP=1


g = natsort.versorted(glob.glob('../../ps_edge_states/B135e-2_F145e-3_2s/B*'))
g = natsort.versorted(glob.glob('../../ps_tmp/ps/B135e-2*'))
g = natsort.versorted(glob.glob('../../ps_LC/B136e-2*'))

tsf = '../../ps_LC/ts_B136e-2_N2e4_Pr1e0_F74e-3_m72_tr1e3_00019_ps'
d   = memmap(tsf,offset=256,dtype=double)
t   = d[::18]
dt  = t[1]-t[0]
del d,t

BUOY     = 2e4

TS_P_SAMP= SKIP*10

om       = 1.35
om       = 1.36

out_base = 'o135e-2_a7431764629e-11'
out_base = 'o136e-2_a74e-3'

print(len(g))

LEFT_EXTREME_CMAP    = myBlues
LEFT_INTERIOR_CMAP   = myBlWh
RIGHT_INTERIOR_CMAP  = myWhRd
RIGHT_EXTREME_CMAP   = myReds

# ======================================================================
# SUBROUTINES ----------------------------------------------------------
def check_file(rec):
  if not os.path.exists(rec):
    print('Reset data file {:s} not found, quitting')
    sys.exit(1)
  return rec

def read_vel(fheader,pdt,udt,pcount=2):
  vel = fromfile(fheader,dtype=udt,count=1)
  pad = fromfile(fheader,dtype=pdt,count=pcount)
  return vel[0].astype(double).T

def reader(f):
  hdt = dtype('(5)i4, (7)f8, i4, f8, (2)i4') # header data type
  pdt = dtype('i4') # padding data type
  asp = 2.          # computational and physical domain aspect ratio
  with open(f,'rb') as fh:
    header = fromfile(fh,dtype=hdt,count=1)
    M,N = header[0][0][1:3]
    t   = header[0][3]
    omega = header[0][1][4]
    buoy_f= header[0][1][1]
    udt = dtype('({:d},{:d}) f8'.format(M+1,N+1))
    uold= read_vel(fh,pdt,udt)
    u   = read_vel(fh,pdt,udt)
    wold= read_vel(fh,pdt,udt)
    w   = read_vel(fh,pdt,udt)
    Told= read_vel(fh,pdt,udt)
    T   = read_vel(fh,pdt,udt,pcount=1)
  Dx,x = cheb(M); x /= asp; Dx *= asp
# Dz,z = cheb(N); z /= asp; Dz *= asp; Dz = Dz.T

# # Compute Additional fields
# u_x,w_x  = Dx@u,Dx@w                   # Deformation ...
# u_z,w_z  = u@Dz,w@Dz                   # ... Tensor
# u_xx,w_xx= Dx@u_x,Dx@w_x               # Second order vel ...
# u_zz,w_zz= u_z@Dz,w_z@Dz               # ... Derivatives
  T_x  = Dx@T                            # Temp. horz. grad.
# T_z  = T@Dz                            # Temp. vert. grad.
# Th   = T - outer(ones(len(z)),z)       # Perturbation Temperature
# Th_z = Th@Dz                           # Pert. Temp. vert. grad. (T_z - 1)
# eta  = u_z - w_x                       # Vorticity
# JDet = u_x * w_z - w_x * u_z           # Deformation Determinant
# uLap = u * (u_xx+u_zz) + w*(w_xx+w_zz) # \vec{u} Lap \vec{u}
# uTh  = u * T_x + w * Th_z              # pert vel cdot grad Th
# X,Z  = meshgrid(x,z,indexing='ij')
# aspect = x.max()/z.max()
  retd = {
#   'x'  : x,
#   'z'  : z,
#   'X'  : X,
#   'Z'  : Z,
#   'u'  : u,
#   'w'  : w,
#   'T'  : T,
#   'Th' : Th,
    'Tx' : T_x,
#   'Tz' : T_z,
#   'Thz': Th_z,
#   'JD' : JDet,
#   'uL' : uLap,
#   'uTh': uTh,
#   'eta': eta,
#   'aspect' : aspect,
#   'x_clip' : 0, #int(       M*CLIPPING*1e-2),
#   'z_clip' : 0, #int(aspect*N*CLIPPING*1e-2),
#   't'             : t,
#   'omega'         : omega,
#   'buoyancy_freq' : buoy_f,
  }
  return retd

def main(f):
  r = reader(f)
  print('READ -- {:s}'.format(f))
  Tx= r['Tx']
  return Tx

def get_POD(Q,K=20):
  N     = len(Q[0])
  Dx,x  = cheb(N-1)
  X,Z   = meshgrid(x,x,indexing='ij')
  L     = array([ q.flatten() for q in Q]).T
  L1    = L[:, :-1]
  L2    = L[:,1:  ]
  print(r'L created')
  U,S,VH= svd(L1)
  u,s,vh= U[:,:K],S[:K],VH[:K] # Truncate
  print(r'SVD computed and truncated to rank {:d}'.format(K))
  w     = vh @ L1.T # rows are modes
  print(r'PCA computed')
  st    = (u.T @ L2 @ vh.T @ diag(1./s))
  print(r'\tilde{S} computed')
  l,p   = eig(st)
  print(r'eigen \tilde{S} computed')
  b     = (u @ p).T # rows are modes
  return S,l,X,Z,w,b

def plot_mode(X_full,Z_full,Mode_full,figname='out',R=0):
  slice_ind = si = arange(R,len(X_full)-R)
  X,Z,U     = X_full[si][si],Z_full[si][si],Mode_full[si][si]
  imi,ima   = U[10:-10,10:-10].min(),U[10:-10,10:-10].max()
  gmi,gma   = U.min(),U.max()
  imf = abs(imi) if abs(imi) >= abs(ima) else abs(ima)
  gmf = abs(gmi) if abs(gmi) >= abs(gma) else abs(gma)
  imi,ima   = -imf,imf
  gmi,gma   = -gmf,gmf
  fig = figure(1,figsize=(8,8)); clf()
  ax  = Axes(fig,[0.,0.,1.,1.])
  ax.set_axis_off()
  fig.add_axes(ax)
  print(figname,gmi,imi,ima,gma)
  if gmi < imi:
    contourf(X,Z,U,levels=linspace(gmi,imi,3),cmap=LEFT_EXTREME_CMAP)
    contour (X,Z,U,levels=linspace(gmi,imi,3),colors='k',linestyles='-')
  if ima < gma:
    contourf(X,Z,U,levels=linspace(ima,gma,3),cmap=RIGHT_EXTREME_CMAP)
    contour (X,Z,U,levels=linspace(ima,gma,3),colors='k',linestyles='-')
  contourf(X,Z,U,levels=linspace(imi,0,7),cmap=LEFT_INTERIOR_CMAP)
  contour (X,Z,U,levels=linspace(imi,0,7),colors='k',linestyles='-')
  contourf(X,Z,U,levels=linspace(0,ima,7),cmap=RIGHT_INTERIOR_CMAP)
  contour (X,Z,U,levels=linspace(0,ima,7),colors='k',linestyles='-')
  if figname is not None:
    savefig('fig/{:s}.png'.format(figname))
  return None

mkl_set_num_threads(1)
h = g[::SKIP]
Pool = mp.Pool(processes=NPROC)
R    = Pool.map(main,h)
Pool.close()
print(len(h))

mkl_set_num_threads(28)
S,l,X,Z,w,b = get_POD(R)

print('Saving PCA Spectrum')
S_out = c_[ 1+arange(len(S)), S/S[0], S/norm(S) ]
NF    = 2
savetxt('{:s}_POD_sigma{:d}.txt'.format(out_base,len(h)-1),
  S_out,
  comments='',
  fmt='%6d'+ NF*' %10.5e',
  header=('{:6s}' + NF*' {:10s}').format('mode','sigma','energy'),
)

print('Saving DMD Spectrum')
z     = log(l)/(.5*dt/pi*BUOY*om*TS_P_SAMP)
print(dt)
l_out = c_[ 1+arange(len(l)), l.real, l.imag, abs(l), z.real, z.imag, abs(z) ]
NF    = 6
savetxt('{:s}_DMD_lambda{:d}.txt'.format(out_base,len(h)-1),
  l_out,
  comments='',
  fmt='%6d'+ NF*' %+10.5e',
  header=('{:>6s}' + NF*' {:>12s}').format('mode','real','imag','abs','lreal','limag','labs'),
)

def ppod(k):
    mode = w[k].reshape(X.shape)
    plot_mode(X,Z,mode,figname='{:s}_pod_mode{:03d}'.format(out_base,k))
    print('done pod {:d}'.format(k))
    return None

def pdmd(k):
    mode = b[k].reshape(X.shape)
    plot_mode(X,Z,mode,figname='{:s}_dmd_mode{:03d}'.format(out_base,k))
    print('done dmd {:d}'.format(k))
    return None


print('Plotting PCA Modes')
Pool = mp.Pool(processes=10)
Pool.map(ppod,range(20))
Pool.map(pdmd,range(20))
Pool.close()
 
print(dt)

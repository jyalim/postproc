import sys,os
from hellaPy import *
from cheb import *
from numpy import *
from pylab import *
from glob import glob
import matplotlib.patches as patches
import multiprocessing as mp

mkl_set_num_threads(2)
NPROCS=16

SKIP=0
CLIPPING=10

CONTOUR_OPTIMIZING = True
CONTOUR_OPTIMIZING = False

PROBE_NORM = sys.argv[4]
PROBE_NORM = ''

TIME_NORM  = False
TIME_NORM  = True
#NORMALIZE_SCALAR  = float(sys.argv[2]) if TIME_NORM else 1.
#INT_NORM          = float(sys.argv[3]) # if TIME_NORM 
IMA = float(sys.argv[2])
GMA = float(sys.argv[3])
if TIME_NORM:
  print('NORMALIZING FIELDS BASED ON PRESCRIBED GLOBAL VALUE')
  print('-- interior post norm: {:f}\n-- norm value: {:f}'.format(1,IMA))

PARTITION_MAPPING = True

PHASE_PLOT  = False

PLOT_FIELDS = ['u','w','eta','T','Tx','Tz','Th','Thz','JD','uL','uTh','ReSu','ReSw','ReCu']
SYMLOG_PLOTS= {
  'u'   : False,
  'w'   : False,
  'eta' : False,
  'T'   : False,
  'Tx'  : False,
  'Tz'  : False,
  'Th'  : False,
  'Thz' : False,
  'JD'  : False,
  'uL'  : False,
  'uTh' : False,
  'ReSu': False,
  'ReSw': False,
  'ReCu': False,
}
LINEAR_PLOTS= {
  'u'   : False,
  'w'   : False,
  'eta' : False,
  'T'   : False,
  'Tx'  : False,
  'Tz'  : False,
  'Th'  : False,
  'Thz' : False,
  'JD'  : False,
  'uL'  : False,
  'uTh' : False,
  'ReSu': False,
  'ReSw': False,
  'ReCu': False,
}

REQ_FIELD = sys.argv[4]
LINEAR_PLOTS[REQ_FIELD]=True

NUM_CONTOURS  = 19
NUM_CONTOURS  =  0

ZERO_CONTOUR  = True

NUM_FILLED_CONTOURS = 256
NUM_FILLED_CONTOURS = 90 
NUM_FILLED_CONTOURS = 15

FIG_BASE      = FBASE= fbase = sys.argv[1]
FIG_DIRECTORY = FDIR = fdir  = 'fig' + '/' + FBASE
OUT_FILE_TYPE = OREC = orec  = 'pdf'
OUT_FILE_TYPE = OREC = orec  = 'png'
fdir = fdir.strip()
if '/' == fdir[-1]:
  fdir = fdir[:-1]
print( 'fdir: {:s}'.format(fdir) )
print('sysr:', sys.argv[5:])

# Full Colormap                # Paint Domain | Paint Range
TOTAL_CMAP           = mycm19  #    [-a, a]   | dark blue to dark red
TOTAL_NLVL           = 19      #    19 colors defined
# Partition Colormaps, for b > a and c < d : a,b,c,d > 0
#                              # Paint Domain | Paint Range
#                              # -------------|------------------
LEFT_EXTREME_CMAP    = myBlues #    [-b,-a]   | dark blue, blue
LEFT_INTERIOR_CMAP   = myBlWh  #    [-a, 0]   | blue,cyan,white
INTERIOR_CMAP        = mycm15  #    [-a, a]   | blue,white,red
RIGHT_INTERIOR_CMAP  = myWhRd  #    [ 0, c]   | white,yellow,red
RIGHT_EXTREME_CMAP   = myReds  #    [ c, d]   | red, dark red

LEFT_EXTREME_NLVL    =  3      #    [-b,-a]   | 3 colors defined
LEFT_INTERIOR_NLVL   =  8      #    [-a, 0]   | 8 colors defined
INTERIOR_NLVL        = 15      #    [-a, a]   |15 colors defined
RIGHT_INTERIOR_NLVL  =  8      #    [ 0, c]   | 8 colors defined
RIGHT_EXTREME_NLVL   =  3      #    [ c, d]   | 3 colors defined

TOLERANCE = 1e-8
# ======================================================================
# SUBROUTINES ----------------------------------------------------------
def check_file(rec):
  if not os.path.exists(rec):
    print(f'Reset data file {rec:s} not found, quitting')
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
  Dz,z = cheb(N); z /= asp; Dz *= asp; Dz = Dz.T

  # Compute Additional fields
  u_x,w_x  = Dx@u,Dx@w                   # Deformation ...
  u_z,w_z  = u@Dz,w@Dz                   # ... Tensor
  u_xx,w_xx= Dx@u_x,Dx@w_x               # Second order vel ...
  u_zz,w_zz= u_z@Dz,w_z@Dz               # ... Derivatives
  T_x  = Dx@T                            # Temp. horz. grad.
  T_z  = T@Dz                            # Temp. vert. grad.
  Th   = T - outer(ones(len(z)),z)       # Perturbation Temperature
  Th_z = Th@Dz                           # Pert. Temp. vert. grad. (T_z - 1)
  eta  = u_z - w_x                       # Vorticity
  ReSu = u @ u_x + w @ u_z
  ReSw = u @ w_x + w @ w_z
  ReCu = ReSu@Dz - Dx@ReSw
  JDet = u_x * w_z - w_x * u_z           # Deformation Determinant
  uLap = u * (u_xx+u_zz) + w*(w_xx+w_zz) # \vec{u} Lap \vec{u}
  uTh  = u * T_x + w * Th_z              # pert vel cdot grad Th
  X,Z  = meshgrid(x,z,indexing='ij')
  aspect = x.max()/z.max()
  retd = {
    'x'   : x,
    'z'   : z,
    'X'   : X,
    'Z'   : Z,
    'u'   : u,
    'w'   : w,
    'T'   : T,
    'Th'  : Th,
    'Tx'  : T_x,
    'Tz'  : T_z,
    'Thz' : Th_z,
    'JD'  : JDet,
    'uL'  : uLap,
    'uTh' : uTh,
    'eta' : eta,
    'ReSu': ReSu,
    'ReSw': ReSw,
    'ReCu': ReCu,
    'aspect' : aspect,
    't'             : t,
    'omega'         : omega,
    'buoyancy_freq' : buoy_f,
  }
  if len(PROBE_NORM)>0:
    q = PROBE_NORM
    Q = retd[q]
    gma = abs(Q).max()
    ind = abs(x - .75*.5).argmin()
    ima = abs(Q[ind:-ind,ind:-ind]).max()
    probe_res = ('{:s} {:s}' + 3*' {:+21.7e}').format(q,fbase,ima,gma,gma/ima)
    print(probe_res,file=sys.stderr)
  return retd

def symlognorm(u):
  v = sign(u) * log(1+abs(u))
  return v / abs(v).max()

def header_print(num_files):
  print(72*'=')
  print('NUM FILES: {:d}'.format(num_files))
  print(72*'-')
  if NPROCS > 1:
    print('PARALLEL MODE')
    print('NPROCS: {:d}'.format(NPROCS))
    print('SERIAL MODE')
  print('PLOTTING FIELDS: ' + (len(PLOT_FIELDS)*'{:s} ').format(*PLOT_FIELDS))
# print('CLIPPING:        {:g}%'.format(CLIPPING))
  print(72*'-')
  print('{:^28s} {:^10s} {:^10s} {:^10s} {:^10s}'.format(
      'file (under {:s}/)'.format(fdir),'imi','ima','gmi','gma'
    )
  )
  return None

def get_figname(f,field,label=''):
  return '{:s}/{:s}_{:s}{:s}.{:s}'.format(fdir,f.split('/')[-1],field,label,OUT_FILE_TYPE)

def check_to_plot(f):
  req = [ False ]
  for field in PLOT_FIELDS:
    out_fig = get_figname(f,field)
    do_plot = False
    if LINEAR_PLOTS[field]:
      if not os.path.exists(out_fig):
        do_plot = True
    elif SYMLOG_PLOTS[field]:
      symlog_out_fig = get_figname(f,field,label='_symlog')
      if not os.path.exists(symlog_out_fig):
        do_plot = True
    req.append(do_plot)
  parse_data = True if any(req) else False
  return parse_data

NEG_EXT_CM = myBlues
INT_CM     = mycm15
POS_EXT_CM = myReds
ALL_CM     = mycm19

def mycf(X,Z,Q,out_fig,ima=1,gma=400,fn=10,ax=0.):
  f,a = no_ax_fax(k=fn,fs_base=6)
  mycontourf(X,Z,Q,levels=linspace(-gma,-ima,3 ),cmap=NEG_EXT_CM)
  mycontourf(X,Z,Q,levels=linspace(-ima, ima,15),cmap=INT_CM    )
  mycontourf(X,Z,Q,levels=linspace( ima, gma,3 ),cmap=POS_EXT_CM)
  if Q.min() < 0 and Q.max() > 0:
    contour(X,Z,Q,levels=[0],linestyles='-',colors='#777777')
  savefig(out_fig)
  return None

Dp,xp  = cheb(512); xp = .5*xp.flatten(); Dp = 2*Dp;
Xp,Zp  = meshgrid(xp,xp,indexing='ij')

def main(f):
  data  = reader(f)
  if len(PROBE_NORM):
    return None
  else:
    x,z       = [ data[key] for key in ['x','z'] ]
    plt_field = sys.argv[4]
    int_field = cheb_interp(2*x,2*z,2*xp,2*xp,{plt_field:data[plt_field]})[plt_field]
    mycf(Xp,Zp,int_field/IMA,get_figname(f,plt_field),ima=1,gma=GMA)
    return int_field
  return None


def main_mean(f_dE):
  f,dE = f_dE
  out_fig=get_figname(f,'eta_pert')
  mycf(Xc,Zc,dE/4500,out_fig,ima=1)
  return None

if __name__ == '__main__':
  drecs = []

  for arg in sys.argv[5:]:
    if '*' in arg:
      drec_list = glob(arg)
      for drec in drec_list:
        drecs.append(check_file(drec))
    else:
      drecs.append(check_file(arg))

  frecs = []
  for drec in drecs:
    plt_f = check_to_plot(drec)
    if not plt_f:
      print( '{:28s} ALREADY PLOTTED FIELDS'.format(get_figname(drec,'$').replace('_$','')) )
    else:
      frecs.append(drec)

  drecs = frecs
  drecs.sort()
   
  if not os.path.exists(fdir):
    os.makedirs(fdir)

  header_print(len(drecs))
  if not CONTOUR_OPTIMIZING:
    if NPROCS > 1:
      qq   = sys.argv[4]
      pool = mp.Pool(processes=NPROCS)
      pool.map(main,drecs[SKIP:])
      if len(PROBE_NORM) > 0:
        pool.map(main,drecs[SKIP:])
      else:
        D = pool.map(main,drecs[SKIP:])
        mkl_set_num_threads(NPROCS)
        D = array(D)
        print(D.shape)
        E = mean(D,axis=0)
        ima = abs(E[125:-125,125:-125]).max()
        mycf(Xp,Zp,E/ima,get_figname(drecs[0].replace('_000000',''),qq+'_mean'),ima=1,gma=GMA)
    # save('D_o143e-2_M2N2_Kxz',D,allow_pickle=True)
    # E   = mean(D,axis=0)
    # ETA = cheb_interp(x,x,q,q,{'eta':E})['eta']
    # ima = abs(ETA[125:-125,125:-125]).max()
    # mycf(Xp,Zp,ETA/ima,'mean_flow_o91e-2_M2N1.png')
    # #dE = load('E.npy')
    # dE  = D-E
    # Z  = [ z for z in zip(drecs,dE) ]
    # pool.map(main_mean,Z)
    else:
      for drec in drecs:
        main(drec)
  else:
    main(drecs[24])

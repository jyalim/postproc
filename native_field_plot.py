import sys,os
from hellaPy import *
from cheb import *
from numpy import *
from pylab import *
from glob import glob
import matplotlib.patches as patches
import multiprocessing as mp

mkl_set_num_threads(2)
NPROCS=40


SKIP=0
CLIPPING=10

INTERP_MESH = 1001
INTERP_MESH = 101
INTERP_MESH = 881
FIG_LENGTH = FL = 10 # inches
FIG_LENGTH = FL = 1  # inches
FIG_LENGTH = FL = 8.8 # inches

PROBE_TX_NORM = True
PROBE_TX_NORM = False

TIME_NORM  = False
TIME_NORM  = True
NORMALIZE_SCALAR  = float(sys.argv[2]) if TIME_NORM else 1.
INT_NORM          = float(sys.argv[3]) # if TIME_NORM 
if TIME_NORM:
  print('NORMALIZING FIELDS BASED ON PRESCRIBED GLOBAL VALUE')
  print('-- interior post norm: {:f}\n-- norm value: {:f}'.format(INT_NORM,NORMALIZE_SCALAR))

PARTITION_MAPPING = True

PHASE_PLOT  = False

PLOT_FIELDS = ['u','w','eta','T','Tx','Tz','Th','Thz','JD','uL','uTh']
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
}
LINEAR_PLOTS= {
  'u'   : False,
  'w'   : False,
  'eta' : False,
  'T'   : False,
  'Tx'  :  True,
  'Tz'  : False,
  'Th'  : False,
  'Thz' : False,
  'JD'  : False,
  'uL'  : False,
  'uTh' : False,
}

NUM_CONTOURS  = 19
NUM_CONTOURS  =  0

NUM_FILLED_CONTOURS = 256
NUM_FILLED_CONTOURS = 90 
NUM_FILLED_CONTOURS = 19

FIG_DIRECTORY = FDIR = fdir  = 'fig' + '/' + sys.argv[1]
OUT_FILE_TYPE = OREC = orec  = 'pdf'
OUT_FILE_TYPE = OREC = orec  = 'png'
fdir = fdir.strip()
if '/' == fdir[-1]:
  fdir = fdir[:-1]
print( 'fdir: {:s}'.format(fdir) )
#print('sysr:', sys.argv[1:])

# Full Colormap                # Paint Domain | Paint Range
TOTAL_CMAP           = mycm19  #    [ a, b]   | dark blue to dark red
# Partition Colormaps, for b > a and c < d : a,b,c,d > 0
#                              # Paint Domain | Paint Range
#                              # -------------|------------------
LEFT_EXTREME_CMAP    = myBlues #    [-b,-a]   | dark blue, blue
LEFT_INTERIOR_CMAP   = myBlWh  #    [-a, 0]   | blue,cyan,white
RIGHT_INTERIOR_CMAP  = myWhRd  #    [ 0, c]   | white,yellow,red
RIGHT_EXTREME_CMAP   = myReds  #    [ c, d]   | red, dark red

LEFT_EXTREME_NLVL    =  3      #    [-b,-a]   | 3 colors defined
LEFT_INTERIOR_NLVL   =  8      #    [-a, 0]   | 8 colors defined
RIGHT_INTERIOR_NLVL  =  8      #    [ 0, c]   | 8 colors defined
RIGHT_EXTREME_NLVL   =  3      #    [ c, d]   | 3 colors defined

TOLERANCE = 1e-8
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
  JDet = u_x * w_z - w_x * u_z           # Deformation Determinant
  uLap = u * (u_xx+u_zz) + w*(w_xx+w_zz) # \vec{u} Lap \vec{u}
  uTh  = u * T_x + w * Th_z              # pert vel cdot grad Th
  X,Z  = meshgrid(x,z,indexing='ij')
  aspect = x.max()/z.max()
  if PROBE_TX_NORM:
    probe_res = 'Tx abs max: {:+21.7e}'.format(abs(T_x).max())
    print(probe_res,file=sys.stderr)

  retd = {
    'x'  : x,
    'z'  : z,
    'X'  : X,
    'Z'  : Z,
    'u'  : u,
    'w'  : w,
    'T'  : T,
    'Th' : Th,
    'Tx' : T_x,
    'Tz' : T_z,
    'Thz': Th_z,
    'JD' : JDet,
    'uL' : uLap,
    'uTh': uTh,
    'eta': eta,
    'aspect' : aspect,
    'x_clip' : int(       M*CLIPPING*1e-2),
    'z_clip' : int(aspect*N*CLIPPING*1e-2),
    't'             : t,
    'omega'         : omega,
    'buoyancy_freq' : buoy_f,
  }
  return retd

def symlognorm(u):
  v = sign(u) * log(1+abs(u))
  return v / abs(v).max()

def mycontour(X,Z,V,mi,ma,num_levels):
  lvls = linspace(mi,ma,num_levels)
  contour(X,Z,V,levels=lvls,colors='k',linestyles='solid',linewidths=1,alpha=.7)
  return None

def mycontourf(X,Z,V,mi,ma,cmap,num_levels):
  cc = contourf(X,Z,V,levels=linspace(mi,ma,NUM_FILLED_CONTOURS),cmap=cmap)
  for c in cc.collections:
    c.set_edgecolor('face')
  if NUM_CONTOURS > 0:
    mycontour(X,Z,V,mi,ma,num_levels)
  return None

def plot_loop(X,Z,V,figsize,imi,ima,gmi,gma,params,out):
  fig = figure(1,figsize=figsize); clf()
  ax  = Axes(fig,[0.,0.,1.,1.])
  unscale  = array([ X.min(),X.max(),Z.min(),Z.max() ])
  rescale  = array([
    1.05 if PHASE_PLOT else 1,
    1,
    1,
    1,
  ])
  axbounds = rescale * unscale
  dscale   = abs(axbounds - unscale)
  ax_w     = axbounds[1] - axbounds[0]
  ax_h     = axbounds[3] - axbounds[2]
  ax.set_axis_off(); ax.axis(axbounds); fig.add_axes(ax)
  if PARTITION_MAPPING == True:
    if imi < 0:
      if gmi < imi:
        mycontourf(X,Z,V,gmi,imi,LEFT_EXTREME_CMAP,LEFT_EXTREME_NLVL)
      mycontourf(X,Z,V,imi,0,LEFT_INTERIOR_CMAP,LEFT_INTERIOR_NLVL)
    if ima > 0:
      mycontourf(X,Z,V,0,ima,RIGHT_INTERIOR_CMAP,RIGHT_INTERIOR_NLVL)
      if ima < gma:
        mycontourf(X,Z,V,ima,gma,RIGHT_EXTREME_CMAP,RIGHT_EXTREME_NLVL)
    if NUM_CONTOURS > 0:
      #contour(X,Z,V,levels=linspace(gmi,gma,NUM_CONTOURS),colors='k',linestyles='solid',linewidths=1,alpha=.7)
      contour(X,Z,V,levels=[0],colors='k',linestyles='solid',linewidths=2,alpha=.3)
  else:
    mycontourf(X,Z,V,gmi,gma,TOTAL_CMAP)
  if PHASE_PLOT == True:
    phase = params['phase']
    xval  = .5 * ( axbounds[0] + X.min() )
    # clean artifacts
    ax.add_patch(
      patches.Rectangle(
        (axbounds[0],axbounds[2]), # (x,z)
                        dscale[0], # width
                             ax_h, # height
                        fill=True, # do background fill
              facecolor='#ffffff', # background fill color
                          alpha=1, # opacity
#                    hatch='....', # hatches
                edgecolor='none', # rectangle edge color
                      linewidth=0, # edge line width
      )
    )
    # Top
    ax.add_patch(
      patches.Rectangle(
        (axbounds[0],axbounds[2]), # (x,z)
                        dscale[0], # width
                       .2351*ax_h, # height
                       fill=False, # do background fill
              facecolor='#ffffff', # background fill color
                         alpha=.5, # opacity
                      hatch='...', # hatches
              edgecolor='#111111', # rectangle edge color
                      linewidth=0, # edge line width
      )
    )
    # Bottom
    ax.add_patch(
      patches.Rectangle(
        (axbounds[0],axbounds[3]), # (x,z)
                        dscale[0], # width
                      -.2351*ax_h, # height
                       fill=False, # do background fill
              facecolor='#ffffff', # background fill color
                         alpha=.5, # opacity
                      hatch='...', # hatches
              edgecolor='#111111', # rectangle edge color
                      linewidth=0, # edge line width
      )
    )
    # Origin
    ax.add_patch(
      patches.Rectangle(
   (axbounds[0],.005*axbounds[3]), # (x,z)
                        dscale[0], # width
                 -.01*axbounds[3], # height
                        fill=True, # do background fill
              facecolor='#111111', # background fill color
                         alpha=.5, # opacity
                        hatch='-', # hatches
              edgecolor='#111111', # rectangle edge color
                      linewidth=0, # edge line width
      )
    )
    #plot(xval, .25      ,'bs',ms=9,clip_on=False,alpha=.2)
    #plot(xval,0.0       ,'ks',ms=9,clip_on=False,alpha=.2)
    plot(xval, .25*phase,'s',mec=[0,0,0,.6],mfc='#0087BC33',ms=9,clip_on=False)#,alpha=.3)
  savefig(out)
  print('{:28s} {:+.3e} {:+.3e} {:+.3e} {:+.3e}'.format(out.split('/')[-1],imi,ima,gmi,gma))
  return None

def get_bounds(field,aspect,clips):
  ic,jc = clips
  field_interior = field[ic:-ic,jc:-jc]
  gmi,gma = field.min(),field.max()
  imi,ima = field_interior.min(),field_interior.max()
  # Symmetrize levels
  imf = abs(imi) if abs(imi) >= abs(ima) else abs(ima)
  gmf = abs(gmi) if abs(gmi) >= abs(gma) else abs(gma)
  if TIME_NORM:
    imf = INT_NORM
    gmf = 1.00
  imi,ima   = -imf,imf
  gmi,gma   = -gmf,gmf
  return imi,ima,gmi,gma

def plot_head(X,Z,U,field_key,aspect,out_fig,clips,params):
  UN    = (U / NORMALIZE_SCALAR) if TIME_NORM else U
  UN[abs(UN)<TOLERANCE] = 0
  imi,ima,gmi,gma = get_bounds(UN,aspect,clips)
  fsize = (FIG_LENGTH,aspect*FIG_LENGTH)
  if LINEAR_PLOTS[field_key]:
    if not os.path.exists(out_fig):
      plot_loop(X,Z,UN,fsize,imi,ima,gmi,gma,params,out_fig)
    else:
      print('{:28s} ALREADY PLOTTED'.format(out_fig))
  V = symlognorm(UN)
  V[abs(V)<TOLERANCE] = 0
  imi,ima,gmi,gma = get_bounds(V,aspect,clips)
  if SYMLOG_PLOTS[field_key]:
    symlog_out_fig = out_fig[:-4]+'_symlog'+out_fig[-4:]
    if not os.path.exists(symlog_out_fig):
      plot_loop(X,Z,V,fsize,imi,ima,gmi,gma,params,symlog_out_fig)
    else:
      print('{:28s} ALREADY PLOTTED'.format(symlog_out_fig))
  return None

def header_print(num_files):
  print(72*'=')
  print('NUM FILES: {:d}'.format(num_files))
  print(72*'-')
  if NPROCS > 1:
    print('PARALLEL MODE')
    print('NPROCS: {:d}'.format(NPROCS))
    print('SERIAL MODE')
  print('PLOTTING FIELDS: ' + (len(PLOT_FIELDS)*'{:s} ').format(*PLOT_FIELDS))
  print('CLIPPING:        {:g}%'.format(CLIPPING))
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

def main(f):
  data  = reader(f)
  asp   = data['aspect']
  x,z   = data['x'],data['z']
  X,Z   = data['X'],data['Z']
  clips = [ data['x_clip'],data['z_clip'] ]
  xp    = linspace(-.5,.5,INTERP_MESH)
  XP,ZP = meshgrid(xp,xp,indexing='ij')
  params= {
    'phase' : cos(data['omega']*data['t']*data['buoyancy_freq']),      
  }
  cheb_fields= { field : data[field] for field in PLOT_FIELDS }
  int_fields = cheb_interp(2*x,2*z,2*xp,2*xp,cheb_fields)
  for field in PLOT_FIELDS:
    plot_head(XP,ZP,int_fields[field],field,asp,get_figname(f,field),clips,params)
  return None

if __name__ == '__main__':
  drecs = []

  for arg in sys.argv[4:]:
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

  if NPROCS > 1:
    pool = mp.Pool(processes=NPROCS)
    pool.map(main,drecs[SKIP:])
  else:
    for drec in drecs:
      main(drec)

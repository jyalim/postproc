from hellaPy import *
from cheb import *
from numpy import *
from pylab import *
import multiprocessing as mp
import glob,natsort
import sys

NPROCS   = 10
NUM_COLS = 18

def get_basename(f):
  return f.split('/')[-1]

def get_token(f,token):
  bn = get_basename(f)
  return bn.split(token)[1].split('_')[0]

def get_fields(f):
  d = memmap(f,dtype=double,offset=256)
  D = array(d).reshape(len(d)//NUM_COLS,NUM_COLS).T
  del d
  t,E,ET,ETh,u0,u1,u2,v0,v1,v2,T0,T1,T2,Kx,Kz,Kxz,Rpi,EP = D
  basename      = get_basename(f)
  omega_str     = get_token(f, 'B'); omega     = float(omega_str)
  alpha_str     = get_token(f, 'F'); alpha     = float(alpha_str)
  buoy_freq_str = get_token(f, 'N'); buoy_freq = float(buoy_freq_str)
  Pr_str        = get_token(f,'Pr'); Pr        = float(Pr_str)
  Gr  = buoy_freq**2
  # Normalizations
  t  *= omega*buoy_freq/(2*pi)
  E  /= Gr
  ET *= 24
  for velo in u0,u1,u2,v0,v1,v2:
    velo /= buoy_freq
  # Save fields in dictionary
  fields = {
    'f'  : f,
    'bn' : basename,
    'oms': omega_str,
    'als': alpha_str,
    'om' : omega,
    'al' : alpha,
    'Ns' : buoy_freq_str,
    'N'  : buoy_freq,
    'Gr' : Gr,
    'Pr' : Pr,
    't'  : t,      
    'E'  : E,
    'ET' : ET,
    'ETh': ETh,
    'u0' : u0,
    'u1' : u1,
    'u2' : u2,
    'v0' : v0,
    'v1' : v1,
    'v2' : v2,
    'T0' : T0,
    'T1' : T1,
    'T2' : T2,
    'Kx' : Kx,
    'Kz' : Kz,
    'Kxz': Kxz,
    'Rpi': Rpi,
    'EP' : EP,
  }
  print('parsed: {:s}'.format(f))
  return fields

def pltter(t,ET,om_s,al_s,out):
  figure(1,figsize=(10,6)); clf()
  plot(t[::200],ET[::200],'k-')
  title(r'$\omega=${:s},$\alpha=${:s}'.format(om_s,al_s))
  xlabel(r'Number of Forcing Periods')
  ylabel(r'$E_T$')
  savefig(out)
  print('plottd: {:s}'.format(out))
  return None

def main(f):
  fields = get_fields(f)
  t,bn,ET,om_s,al_s = [ fields[key] for key in ['t','bn','ET','oms','als']]
  out = 'fig/{:s}_ET.png'.format(bn)
  pltter(t,ET,om_s,al_s,out)
  return None

G = natsort.realsorted(glob.glob(sys.argv[1]))

p = mp.Pool(processes=NPROCS)
p.map(main,G)

savetxt('var/file_list',G,fmt='%s',comments='')

#for g in G:
#  main(g)

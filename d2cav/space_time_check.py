from hellaPy import *
from numpy import *
from pylab import *
import glob,natsort,sys

def get_data(f):
  d = memmap(f,dtype=double,offset=256,mode='r')
  K = 100*(len(d)//1800)
  J = K-5000*100
  u0= d[4::18][J:K+1][::100]/2e4
  u1= d[5::18][J:K+1][::100]/2e4
  w0= d[7::18][J:K+1][::100]/2e4
  w1=-d[8::18][J:K+1][::100]/2e4
  s = 0 if u0[0] > 0 else 1
  u0= u0[s  ::2]
  w0= w0[s  ::2]
  u1= u1[s^1::2]
  w1= w1[s^1::2]
  return u0,w0,u1,w1

def get_token(f,token):
  return f.split(token)[-1].split('_')[0]

def main(f):
  if '_F' in f:
    om,al = [ get_token(f,token) for token in ['B','F'] ]
  else:
    om,al = [ get_token(f,token) for token in ['_o','_a'] ]
  u0,w0,u1,w1 = get_data(f)
  figure(1,figsize=(8,8)); clf()
  plot(u0,w0,'o',mec='none',mfc='k',ms=1)
  plot(u1,w1,'o',mec='none',mfc='r',ms=1)
  xlabel(r'$u$')
  ylabel(r'$w$')
  title(r'$(\omega,\alpha)=$'+f'({om:s},{al:s})',y=1.10)
  savefig(f'fig/space_time_monitor/o{om:s}_a{al:s}.png')
  print(f'plotted: {f:s}')
  return None

G = natsort.realsorted(glob.glob(sys.argv[1]))
for g in G:
  main(g)

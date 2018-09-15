from hellaPy import *
from numpy import *
from pylab import *
import glob,natsort,sys
import multiprocessing as mp
import pandas as pd
import csv

def get_data(f):
  d = memmap(f,dtype=double,offset=256,mode='r')
  K = 100*(len(d)//1800)
  J = K-4000*100
  u0= d[4::18][J:K+1]/2e4
  u1= d[5::18][J:K+1]/2e4
  w0= d[7::18][J:K+1]/2e4
  w1= d[8::18][J:K+1]/2e4
  return u0.std(),w0.std(),u1.std(),w1.std()

def get_token(f,token):
  return f.split(token)[-1].split('_')[0]

def main(f):
  if '_F' in f:
    om,al = [ get_token(f,token) for token in ['B','F'] ]
  else:
    om,al = [ get_token(f,token) for token in ['_o','_a'] ]
  u0,w0,u1,w1 = get_data(f)
  print(f'PARSED: {f:s}')
  return [om,al,u0,w0,u1,w1]

F  = sys.argv[1]
out= sys.argv[2]
G  = natsort.realsorted(glob.glob(F))
p  = mp.Pool(2)
D  = p.map(main,G)
df = pd.DataFrame(D,columns=['omega','alpha','u0','w0','u1','w1'])
ind= argsort(df.alpha.astype(float))
df.iloc[ind].to_csv(out,sep=' ',float_format='%16.7e',index_label=False,index=False,quoting=csv.QUOTE_NONE,escapechar=' ')

#for g in G:
#  D=main(g)

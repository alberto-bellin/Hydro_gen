#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 09:46:18 2018

@author: alberto
"""

f = open('search_volume.txt', 'w')

k=0
jfin=4
ijcount=jfin*3
n1=ijcount*2+ijcount-2
nn2=4
nn1=3


k=0
for j in range(1,n1+1):
    k2=int((j-1)/ijcount)
    jp=j-k2*ijcount
    i2=int((jp-1)/jfin)
    for i in range(1,j+1):
        k1=int((i-1)/ijcount)
        ip=i-k1*ijcount
        i1=int((ip-1)/jfin)
        ilag=abs(i2-i1)+1
        jlag=abs(jp-ip-(i2-i1)*jfin)+1
        klag=k2-k1+1
        k=k+1
#        s='k='+repr(k)+' ilag='+repr(ilag)+' jlag='+repr(jlag)+' klag='+repr(klag)
        f.write('k='+str(k)+' ilag='+str(ilag)+' jlag='+str(jlag)+' klag='+str(klag)+'\n')
        
f.write('-------------'+'\n')


#f.close()

#  Known term

#f = open('search_volume1.txt', 'w')

k=0
k2=int(n1/ijcount)
jp=n1+1-k2*ijcount
i2=int((jp-1)/jfin)    
for i in range(1,n1+1):
    k1=int((i-1)/ijcount)
    ip=i-k1*ijcount
    i1=int((ip-1)/jfin)
    ilag=abs(i2-i1)+1
    jlag=abs(jp-ip-(i2-i1)*jfin)+1
    klag=k2-k1+1
    k=k+1
    str(i) 
    f.write('k='+str(k)+' ilag='+str(ilag)+' jlag='+str(jlag)+' klag='+str(klag)+'\n')
f.close()    
    

#!/usr/bin/env python
#-*- coding:utf-8 -*-
from __future__ import print_function, division
(T,F)=(True,False)

sw=0 #0: ham_r, 1: .input, 2:_hr.dat
name='000AsP.input'
brav='P' #set bravais lattice Initial
so=T #switch g-vector (off diagonal part of so hamiltonian)
eps=1.0e-6

import numpy as np
def read_ham(sw):
    def import_hop():
        rvec=np.loadtxt('irvec.txt')
        nr=rvec[:,0].size
        tmp=np.array([complex(float(tp[0]),float(tp[1])) for tp in [f.strip(' ()\n').split(',') for f in open('ham_r.txt','r')]])
        no=int(np.sqrt(tmp.size/nr))
        ham_r=tmp.reshape(nr,no,no)
        ndegen=(np.loadtxt('ndegen.txt') if True else np.ones(nr))
        return(rvec,ndegen,ham_r,no,nr)
    def import_out(fname):
        data=np.loadtxt(fname)
        con=(data[:,:3]==data[0,:3]).prod(axis=1).sum()
        no,nr =int(np.sqrt(con)),data[:,0].size//con
        rvec=data[:nr,:3]
        ham_r=(data[:,3]+1j*data[:,4]).reshape(nr,no,no)
        ndegen=np.ones(nr)
        return(rvec,ndegen,ham_r,no,nr)
    def import_hr(name):
        tmp=[f.split() for f in open('%s_hr.dat'%name,'r')]
        no, nr=int(tmp[1][0]), int(tmp[2][0])
        c2,tmp1=3,[]
        while not len(tmp1)==nr:
            tmp1.extend(tmp[c2])
            c2=c2+1
        ndegen=np.array([int(t) for t in tmp1])
        tmp1=[[float(t) for t in tp] for tp in tmp[c2:]]
        tmp=np.array([complex(tp[5],tp[6]) for tp in tmp1])
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        ham_r=tmp.reshape(nr,no,no)
        return(rvec,ndegen,ham_r,no,nr)
    def import_Hopping():
        tmp=[f.split() for f in open('Hopping.dat','r')]
        axis=np.array([[float(tp) for tp in tpp] for tpp in tmp[1:4]])
        no,nr=int(tmp[4][0]),int(tmp[4][1])
        ndegen=np.array([1]*nr)
        tmp1=[[float(t) for t in tp] for tp in tmp[7+no:]]
        rvec=np.array([tmp1[no*no*i][:3] for i in range(nr)])
        tmp=np.array([complex(tp[8],tp[9]) for tp in tmp1])
        ham_r=tmp.reshape(nr,no,no)
        return(rvec,ndegen,ham_r,no,nr)

    return(import_hop() if sw==0 else import_out(name) 
           if sw==1 else import_hr(name) if sw==2 else import_Hopping())

def check_ham(ham_r,rvec,f,sw=True):
    count=0
    for r1,hm1 in zip(rvec,ham_r):
        if(abs(hm1).sum()>eps):
            r=(-r1 if sw else r1)
            for r2,hm2 in zip(rvec,ham_r):
                if(abs(r-r2).sum()<1.0e-3):
                    tmp=f(hm1,hm2).sum()/hm1.size
                    if(tmp>eps):
                        count=count+1
                    break
    return(count)

def check_hermite(ham_r,rvec):
    f=lambda a,b:abs(a-b.T.conjugate()) #ck t(r)=t^+(-r)
    count=check_ham(ham_r,rvec,f)
    print(('' if(count==0) else 'not ')+'Hermite'+('' if(count==0) else '\n%d'%count))

def check_SRS(ham_r,rvec): #proper only single site
    f=lambda a,b:abs(a-b) #ck t(r)=t(-r)
    count=check_ham(ham_r,rvec,f)
    print('SRS'+('' if(count==0)else ' breaking\n%d'%count))

def check_TRS(ham_r,rvec):
    sw_f=False
    if(so):
        '''
        K=(Huu+Hdd)/2,gx=(Hud+Hdu)/2,gy=i(Hud-Hdu)/2,gz=(Huu-Hdd)/2
        if Hamiltonian has TRS K(r)=K*(r) and g_i(r)=-g*_i(r) (i=x,y,z)
        '''
        no=int(len(ham_r[0])/2)
        f=lambda a,b:abs(a-b.conjugate())
        hm=(ham_r[:,:no,:no]+ham_r[:,no:2*no,no:2*no])*0.5  #K
        count=check_ham(hm,rvec,f,sw_f) #ck K(r)=K*(r)
        f=lambda a,b:abs(a+b.conjugate())
        hm=(ham_r[:,:no,no:2*no]+ham_r[:,no:2*no,:no])*0.5  #gx
        count=count+check_ham(hm,rvec,f,sw_f) #ck gx(r)=-gx*(r)
        hm=(ham_r[:,:no,no:2*no]-ham_r[:,no:2*no,:no])*0.5j #gy
        count=count+check_ham(hm,rvec,f,sw_f) #ck gy(r)=-gy*(r)
        hm=(ham_r[:,:no,:no]-ham_r[:,no:2*no,no:2*no])*0.5  #gz
        count=count+check_ham(hm,rvec,f,sw_f) #ck gz(r)=-gz*(r)
    else:
        #f=lambda a,b:abs(a-b.T)
        f=lambda a,b:abs(a-b.conjugate())
        count=check_ham(ham_r,rvec,f,sw_f)

    print('TRS'+('' if(count==0)else ' breaking\n%d'%count))

def check_imag(ham_r,rvec):
    count=0
    for r,hm in zip(rvec,ham_r):
        if(abs(hm.imag).sum()/hm.size>eps):
            count=count+1
    print('hopping matrices are '+('real' if(count==0)else 'complex'))
    return (True if(count==0) else False)

def check_onsite_energy(rvec,ham_r,no):
    cvec=np.array([0.0,0.0,0.0])
    print('check on-site energy...')
    for r,hm in zip(rvec,ham_r):
        if(abs(r-cvec).sum()<eps):
            onsite_energy=hm.diagonal()
            break
    for j,oe in enumerate(onsite_energy):
        print('orbital num:%3d energy: %f'%(j+1,oe.real))
    return(onsite_energy)

def find_ham_r(rvec,ham_r,ri):
    for r,ham in zip(rvec,ham_r):
        if(abs(r-ri).sum()<eps):
            hm=ham
            break
    return(hm)

def print_ham_r(hm,flag):
    print('   ',end='')
    for i in range(len(hm)):
        print(('     %2d      '%(i+1) if flag 
               else '%11s %2d %11s'%('',(i+1),'')) ,end='')
    print('')
    for hmm in hm:
        print('%2d '%(i+1),end='')
        for hmmm in hmm:
            if flag:
                print(' %11.4e,'%round(hmmm.real,5),end='')
            else:
                print('(%11.4e,%11.4e),'%(round(hmmm.real,5),round(hmmm.imag,5)),end='')
        print('')
    print('')

def print_gvec(hm,f,no0):
    print('   ',end='')
    for i in range(no0):
        print('    %2d      '%(i+1),end='')
    print('')
    for i in range(no0):
        print('%2d '%(i+1),end='')
        for j in range(no0):
            hmm=f(hm,i,j)
            if abs(hmm.imag)<1.e-6:
                print(' %9.2e, '%round(hmm.real,5),end='')
            elif abs(hmm.real)<1.e-6:
                print('%9.2ej, '%round(hmm.imag,5),end='')
            else:
                print('(%9.2e,%9.2e) '%(round(hmm.real,5),round(hmm.imag,5)),end='')
        print('')
    print('')

def main():
    (rvec,ndegen,ham_r,no,nr)=read_ham(sw)
    rflag=check_imag(ham_r,rvec)
    check_hermite(ham_r,rvec)
    check_SRS(ham_r,rvec)
    check_TRS(ham_r,rvec)
    print('nr = %d no = %d'%(nr,no))
    onsite_energy=check_onsite_energy(rvec,ham_r,no)

    """usually chech_r is [-1.,1.,0.],[-1.,0.,0.],[1.,1.,0.],[0.,-1.,-1.],[0.,0.,-1.],[1.,1.,1.]"""

    check_r=[[1.,0.,0.],[-1.,0.,0.]] #,[1.,1.,0.]]
    for r in check_r:
        hm=find_ham_r(rvec,ham_r,r)
        if so:
            b0=lambda a,x,y:(a[x][y]+a[x+int(no*0.5)][y+int(no*0.5)])*0.5
            print('kinetic term')
            print_gvec(hm,b0,int(no*0.5)) #ek

            b0=lambda a,x,y:(a[x][y+int(no*0.5)]+a[x+int(no*0.5)][y])*0.5
            print('gx')
            print_gvec(hm,b0,int(no*0.5)) #gx

            b0=lambda a,x,y:(a[x][y+int(no*0.5)]-a[x+int(no*0.5)][y])*0.5j
            print('gy')
            print_gvec(hm,b0,int(no*0.5)) #gy

            b0=lambda a,x,y:(a[x][y]-a[x+int(no*0.5)][y+int(no*0.5)])*0.5
            print('gz')
            print_gvec(hm,b0,int(no*0.5)) #gz
        else:
            print_ham_r(hm,rflag)

if __name__=='__main__':
    main()

#==============================================================================#
# This program is hopping symmetry checker script                              #
#                                                                              #
# Copyright (c) 2018-2019  K. Suzuki                                           #
#==============================================================================#
#                          The MIT License (MIT)                               #
#==============================================================================#
#Permission is hereby granted, free of charge, to any person obtaining a copy  #
#of this software and associated documentation files (the "Software"), to deal #
#in the Software without restriction, including without limitation the rights  #
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
#copies of the Software, and to permit persons to whom the Software is         #
#furnished to do so, subject to the following conditions:                      #
#                                                                              #
#The above copyright notice and this permission notice shall be included in all#
#copies or substantial portions of the Software.                               #
#                                                                              #
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE #
#SOFTWARE.                                                                     #
#==============================================================================#

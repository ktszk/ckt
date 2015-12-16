from __future__ import print_function
sw=2
name='FeS'

def import_hop():
    import math
    tmp=[f.split() for f in open('irvec.txt','r')]
    rvec=[[float(tp[0]),float(tp[1]),float(tp[2])] for tp in tmp]
    nr=len(rvec)
    ndegen=([float(f) for f in open('ndegen.txt','r')] if True else [1]*nr)
    no=int(math.sqrt(sum(1 for line in open('ham_r.txt','r'))/nr))
    tmp=[f.strip(' ()\n').split(',') for f in open('ham_r.txt','r')]
    tmp1=[complex(float(tp[0]),float(tp[1])) for tp in tmp]
    ham_r=[[[tmp1[i*no*no+j*no+k] for k in range(no)] for j in range(no)] for i in range(nr)]
    return(rvec,ndegen,ham_r,no,nr)

def import_out(fname):
    import math
    tmp=[f.split() for f in open(fname,'r')]
    tmp1=[[float(tp) for tp in tpp] for tpp in tmp]
    tmp=[complex(tp[3],tp[4]) for tp in tmp1]
    r0=tmp1[0][:3]
    con=sum(1 if rr[:3]==r0 else 0  for rr in tmp1)
    no=int(math.sqrt(con))
    nr=len(tmp1)/con
    rvec=[tp[:3] for tp in tmp1[:nr]]
    ham_r=[[[tmp[k+j*nr+i*nr*no] for k in range(no)] for j in range(no)] for i in range(nr)]
    ndegen=[1]*nr
    return(rvec,ndegen,ham_r,no,nr)

def import_hr(name):
    tmp=[f.split() for f in open('%s_hr.dat'%name)]
    no=int(tmp[1][0])
    nr=int(tmp[2][0])
    c1=0; c2=3
    while not c1==nr:
        c1=c1+len(tmp[c2])
        c2=c2+1
    tmp1=[[int(t) for t in tp] for tp in tmp[3:c2]]
    ndegen=[]
    for tp in tmp1:
        ndegen=ndegen+tp
    tmp1=[[float(t) for t in tp] for tp in tmp[c2:]]
    tmp=[complex(tp[5],tp[6]) for tp in tmp1]
    rvec=[tmp1[no*no*i][:3] for i in range(nr)]
    ham_r=[[[tmp[k+j*no+i*no*no] for k in range(no)] for j in range(no)] for i in range(nr)]
    return(rvec,ndegen,ham_r,no,nr)

def check_ham(ham_r,rvec,f):
    count=0
    for i,r1 in enumerate(rvec):
        r=[-rr for rr in r1]
        for j,r2 in enumerate(rvec):
            if(r==r2):
                for k,hmm in enumerate(ham_r[i]):
                    for l,hm in enumerate(hmm):
                        tmp=f(hm,ham_r[j],k,l)
                        if(abs(hm)>1.0e-6 and tmp>1.0e-8):
                            count=count+1
                            #print(k,l,r)
                            #print(hm,ham_r[j][k][l])
                break
    return(count)

def check_hermite(ham_r,rvec):
    f=lambda a,b,x,y:abs(a-b[y][x].conjugate())/(1 if abs(a)==0 else abs(a))
    count=check_ham(ham_r,rvec,f)
    print(('' if(count==0) else 'not ')+'Hermite')

def check_SRS(ham_r,rvec):
    f=lambda a,b,x,y:abs(a-b[x][y])
    count=check_ham(ham_r,rvec,f)
    print('SRS'+('' if(count==0)else ' breaking\n%d'%count))

def check_periodic(ham_r,rvec):
    f=lambda a,b,x,y:abs(a-b[y][x])
    count=check_ham(ham_r,rvec,f)
    print('Periodic symmetry'+('' if(count==0)else ' breaking\n%d'%count))

def write_ham_r(ham_r,filename):
    f=open(filename,'w')
    for hmm in ham_r:
        for hm in hmm:
            for h in hm:
                f.write(' (%18.15E,%18.15E)\n'%(h.real,h.imag))
    f.close()

def check_energy(rvec,ham_r,no,cvec):
    for i,r in enumerate(rvec):
        if(r==cvec):
            out_energy=[ham_r[i][j][j] for j in range(no)]
            for j,oe in enumerate(out_energy):
                print('orbital num:%3d energy: %f'%(j+1,oe.real))
            return(out_energy)

def check_onsite_energy(rvec,ham_r,no):
    cvec=[0.0,0.0,0.0]
    print('check on-site energy...')
    onsite_energy=check_energy(rvec,ham_r,no,cvec)
    return(onsite_energy)

def check_t1(rvec,ham_r,no):
    cvec=[1.0,0.0,0.0]
    print ('check t1 direct')
    t1=check_energy(rvec,ham_r,cvec)
    return(t1)

def check_t2(rvec,ham_r,no):
    cvec=[1.0,1.0,0.0]
    print ('check t2 direct')
    t2=check_energy(rvec,ham_r,no,cvec)
    return(t2)

def check_inter_energy(rvec,ham_r,no,no2,cvec):
    for i,r in enumerate(rvec):
        if(r==cvec):
            for k in range(no2):
                out_energy=[ham_r[i][j][no+k] for j in range(no)]
                for j,oe in enumerate(out_energy):
                    print('orbital num:%3d to %3d energy: %f'%(j+1,no+k+1,oe.real))
            return(0)

def check_t1_indirect(rvec,ham_r,no):
    cvec=[1.0,0.0,0.0]
    print('check t1 indirect')
    t1=check_inter_energy(rvec,ham_r,10,6,cvec)

def check_t2_indirect(rvec,ham_r,no):
    cvec=[1.0,1.0,0.0]
    print('check t2 indirect')
    t2=check_inter_energy(rvec,ham_r,10,6,cvec)

def mkham(ham_r,no):
    f=open('ham_r3.txt','w')
    orb_list=[1,2,5]
    for hm in ham_r:
        for j in range(no):
            for k in range(no):
                if (j+1 in orb_list) and (k+1 in orb_list):
                    f.write(' (%18.15E,%18.15E)\n'%(hm[j][k].real,hm[j][k].imag))
    f.close()

def mk_non_so_spin_model(ham_r,no):
    f=open('ham_r3.txt','w')
    fc=lambda x,y:x%y+2*int(x/(y*2))
    for hm in ham_r:
        for j in range(2*no):
            for k in range(2*no):
                jj=fc(j,no)
                kk=fc(k,no)
                a=((hm[jj][kk].real,hm[jj][kk].imag) 
                   if (j==jj and k==kk) or (j!=jj and k!=kk) else (0.,0.))
                f.write(' (%18.15E,%18.15E)\n'%a)
    f.close()

def find_ham_r(rvec,ham_r,ri):
    for i,r in enumerate(rvec):
        if(r==ri):
            hm=ham_r[i]
            break
    return(hm)

def print_ham_r(hm,b0):
    for i,hmm in enumerate(hm):
        for j,hmmm in enumerate(hmm):
            if(b0(i) and b0(j)):
                print('(%11.4e,%11.4e) %d %d'%(round(hmmm.real,5),round(hmmm.imag,5),i+1,j+1),end='')
        if(b0(i)):
            print('')
    print('')

def print_gvec(hm,f,no):
    for i in range(no):
        for j in range(no):
            hmm=f(hm,i,j)
            print('(%11.4e,%11.4e) %d %d'%(round(hmm.real,5),round(hmm.imag,5),i+1,j+1),end='')
        print('')
    print('')

def IBSC_unfold(rvec,ham_r):
    inv_list=[3]
    r2=[[r[0]+r[1],r[1]-r[0],r[2]] for r in rvec]
    r3=[[r[0]+r[1],r[1]-r[0]+1,r[2]] for r in rvec]
    rvec2=r2+r3
    ham1=[[h[:5] for h in hm[:5]] for hm in ham_r]
    ham2=[[[ h if (i in inv_list) else -h for i,h in enumerate(hm[5:])] for hm in ham[:5]] for ham in ham_r]
    ham_r2=ham1+ham2
    return(rvec2,ham_r2)

def output_unfold(ham_r,rvec):
    write_ham_r(ham_r,'ham_r2.txt')
    f=open('irvec2.txt','w')
    for r in rvec:
        f.write('%12i %11i %11i\n'%tuple(r))
    f.close()

def restruct_ham_r(ham_r,ndegen):
    for i,(n,hmm) in enumerate(zip(ndegen,ham_r)):
        for j,hm in enumerate(hmm):
            for k,h in enumerate(hm):
                ham_r[i][j][k]=h/n
    return(ham_r)

def main():
    (rvec,ndegen,ham_r,no,nr)=(import_hop() if sw==0 else import_out(name) 
                               if sw==1 else import_hr(name))
    ham_r=restruct_ham_r(ham_r,ndegen)
    #(rvec2,ham_r2)=IBSC_unfold(rvec,ham_r)
    #output_unfold(ham_r2,rvec2)
    #write_ham_r(ham_r,'ham_r2.txt')
    #mkham(ham_r,no)
    #mk_non_so_spin_model(ham_r,no)
    check_hermite(ham_r,rvec)
    check_SRS(ham_r,rvec)
    check_periodic(ham_r,rvec)
    print('nr = %d no = %d'%(nr,no))
    onsite_energy=check_onsite_energy(rvec,ham_r,no)

    #hm=find_ham_r(rvec,ham_r,[0.,0.,0.])
    #b0=lambda x:x < no/2
    #print_ham_r(hm,b0)
    #b0=lambda x:x >= no
    #print_ham_r(hm,b0)

    #hm=find_ham_r(rvec,ham_r,[1.,0.,0.])
    #hm=find_ham_r(rvec,ham_r,[0.,1.,0.])
    #hm=find_ham_r(rvec,ham_r,[0.,1.,1.])
    #b0=lambda x:x < no
    #print_ham_r(hm,b0)

    #b0=lambda a,x,y:(a[x+no/2][y]+a[x][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:-(a[x+no/2][y]-a[x][y+no/2])*0.5j
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:(a[x][y]-a[x+no/2][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:(a[x][y]+a[x+no/2][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

    #hm=find_ham_r(rvec,ham_r,[1.,1.,2.])
    #hm=find_ham_r(rvec,ham_r,[-1.,1.,0.])
    #hm=find_ham_r(rvec,ham_r,[-1.,0.,0.])
    #hm=find_ham_r(rvec,ham_r,[1.,1.,0.])
    #hm=find_ham_r(rvec,ham_r,[0.,-1.,-1.])
    #hm=find_ham_r(rvec,ham_r,[0.,0.,-1.])
    #hm=find_ham_r(rvec,ham_r,[1.,1.,1.])
    #b0=lambda x:x < no
    #print_ham_r(hm,b0)

    #b0=lambda a,x,y:(a[x+no/2][y]+a[x][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:-(a[x+no/2][y]-a[x][y+no/2])*0.5j
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:(a[x][y]-a[x+no/2][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

    #b0=lambda a,x,y:(a[x][y]+a[x+no/2][y+no/2])*0.5
    #print_gvec(hm,b0,no/2)

"""
    onsite_energy_Fe=onsite_energy[0:10]
    onsite_energy_As=onsite_energy[10:]
    aveFe=sum(oe.real for oe in onsite_energy_Fe)*0.1
    aveAs=sum(oe.real for oe in onsite_energy_As)/6.0
    print('diff between Fe and As=%f'%(aveFe-aveAs))

    t=check_t1(rvec,ham_r,no)
    t=check_t2(rvec,ham_r,no)
    check_t1_indirect(rvec,ham_r,no)
    check_t2_indirect(rvec,ham_r,no)
"""
if __name__=='__main__':
    main()

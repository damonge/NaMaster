import numpy as np
import matplotlib.pyplot as plt
import os
import healpy as hp
import pymaster as nmt
import sys

DTOR=np.pi/180

def getmaskapoana(ns,aps,fsk=0.1) :
    #This generates a correctly-apodized mask
    vv=hp.pix2vec(ns,np.arange(hp.nside2npix(ns)));
    cth=vv[0]; th=np.arccos(cth); th0=np.arccos(1-2*fsk); th_apo=aps*DTOR
    print th0/DTOR
    id0=np.where(th>=th0)[0]
    id1=np.where(th<=th0-th_apo)[0]
    idb=np.where((th>th0-th_apo) & (th<th0))[0]
    x=np.sqrt((1-np.cos(th[idb]-th0))/(1-np.cos(th_apo)))
    mask_apo=np.zeros(hp.nside2npix(ns))
    mask_apo[id0]=0.
    mask_apo[id1]=1.
    mask_apo[idb]=x-np.sin(2*np.pi*x)/(2*np.pi)
    return mask_apo

if len(sys.argv)!=9 :
    print "python check_sph_pure.py nside w_cont nsims plotres aposize apotype pureE pureB"
    exit(1)

nside  =  int(sys.argv[1])
w_cont =  int(sys.argv[2])
nsims  =  int(sys.argv[3])
plotres=  int(sys.argv[4])
aposize=float(sys.argv[5])
apotype=      sys.argv[6]
pureE=    int(sys.argv[7])
pureB=    int(sys.argv[8])
ispure_e=False
if pureE!=0 :
    ispure_e=True
ispure_b=False
if pureB!=0 :
    ispure_b=True
if w_cont :
    raise ValueError("Contaminant + purification still not implemented")

#alpha_cont_2=0.1
predir="tests_sph"
os.system("mkdir -p "+predir)
prefix=predir+"/run_pure%d%d_ns%d_cont%d_apo%.2lf"%(pureE,pureB,nside,w_cont,aposize)+apotype
fname_mask=prefix+"_mask"

#This just generates the theory power spectra
data=np.loadtxt('planck1_r0p00_lensedtotCls.dat',unpack=True)

ell=np.arange(3*nside)
cltt=np.zeros(3*nside); clee=np.zeros(3*nside); clbb=np.zeros(3*nside); clte=np.zeros(3*nside); 
cltt[2:]=data[1,:3*nside-2]*2*np.pi/(data[0,:3*nside-2]*(data[0,:3*nside-2]+1.))
clee[2:]=data[2,:3*nside-2]*2*np.pi/(data[0,:3*nside-2]*(data[0,:3*nside-2]+1.))
clbb[2:]=data[3,:3*nside-2]*2*np.pi/(data[0,:3*nside-2]*(data[0,:3*nside-2]+1.))
clte[2:]=data[4,:3*nside-2]*2*np.pi/(data[0,:3*nside-2]*(data[0,:3*nside-2]+1.))
if plotres :
    plt.figure()
    plt.plot(ell,cltt,'r-',label='$TT$')
    plt.plot(ell,clee,'b-',label='$EE$')
    plt.plot(ell,clbb,'g-',label='$BB$')
    plt.plot(ell,clte,'y-',label='$TE$')
    plt.loglog()
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='upper right',frameon=False,fontsize=16,labelspacing=0.1,ncol=2)

if w_cont :
    '''
    tilt_fg=-2.0
    l0_fg=100.
    clttfg=1E-5*((ell+10.)/(l0_fg+10.))**tilt_fg
    cleefg=5E-7*((ell+30.)/(l0_fg+30.))**tilt_fg
    cltefg=0.9*np.sqrt(clttfg*cleefg)
    clbbfg=0.5*cleefg
    clttfg[0]=0; cleefg[0]=0; clbbfg[0]=0; cltefg[0]=0; 
    if plotres :
        plt.plot(ell,alpha_cont_0*alpha_cont_0*clttfg,'r--',label='${\\rm FG},\\,TT$')
        plt.plot(ell,alpha_cont_0*alpha_cont_2*cltefg,'g--',label='${\\rm FG},\\,TE$')
        plt.plot(ell,alpha_cont_2*alpha_cont_2*cleefg,'b--',label='${\\rm FG},\\,EE$')
        plt.plot(ell,alpha_cont_2*alpha_cont_2*clbbfg,'y--',label='${\\rm FG},\\,BB$')
    '''
    raise ValueError("Can't do purification and deprojection yet")


np.random.seed(1001)
fsky=0.2
rholes=1.
if not os.path.isfile(fname_mask+'.fits') :
    print "Generating mask"
    mask=getmaskapoana(nside,aposize,fsk=fsky)
    hp.write_map(fname_mask+".fits",mask)
mask=hp.read_map(fname_mask+".fits")
if plotres :
    hp.mollview(mask)

if w_cont :  #Not ready yet
    if not os.path.isfile(prefix+"_contaminants.fits") :
        fgt,fgq,fgu=nmt.synfast_spherical(nside,[clttfg,cleefg,clbbfg,cltefg],pol=True)
        hp.write_map(prefix+"_contaminants.fits",[fgt,fgq,fgu])
    else :
        fgt,fgq,fgu=hp.read_map(prefix+"_contaminants.fits",field=[0,1,2],verbose=False)

#Binning scheme
d_ell=int(1./fsky)
b=nmt.NmtBin(nside,nlb=d_ell)

#Generate some initial fields
print " - Res: %.3lf arcmin. "%(np.sqrt(4*np.pi*(180*60/np.pi)**2/hp.nside2npix(nside)))
def get_fields() :
    mppt,mppq,mppu=nmt.synfast_spherical(nside,[cltt,clee,clbb,clte],pol=True)
    if w_cont : #Not ready yet
        mppq+=alpha_cont_2*fgq
        mppu+=alpha_cont_2*fgu
        ff2=nmt.NmtField(mask,[mppq,mppu],[[fgq,fgu]])
    else :
        ff2=nmt.NmtField(mask,[mppq,mppu],purify_e=ispure_e,purify_b=ispure_b,n_iter_mask_purify=10)
    return mppq,mppu,ff2
mpq,mpu,f2=get_fields()

if plotres :
    hp.mollview((mpq*mask).flatten(),title='$Q$')
    hp.mollview((mpu*mask).flatten(),title='$U$')

#Compute deprojection bias
if w_cont : #Not ready yet
    clb22=nmt.deprojection_bias(f2,f2,[clee,0*clee,0*clbb,clbb])
else :
    clb22=None;

#Use initial fields to generate coupling matrix
w22=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w22.dat") :
    print "Computing 22"
    w22.compute_coupling_matrix(f2,f2,b)
    w22.write_to(prefix+"_w22.dat");
else :
    w22.read_from(prefix+"_w22.dat")

#Generate theory prediction
cl22_th=w22.decouple_cell(w22.couple_cell([clee,0*clee,0*clbb,clbb]))
np.savetxt(prefix+"_cl_th.txt",
           np.transpose([b.get_effective_ells(),
                         cl22_th[0],cl22_th[1],cl22_th[2],cl22_th[3]]))

#Compute mean and variance over nsims simulations
cl22_all=[]
for i in np.arange(nsims) :
    if i%100==0 :
        print "%d-th sim"%i

    if not os.path.isfile(prefix+"_cl_%04d.txt"%i) :
        mpq,mpu,f2=get_fields()
        cl22=w22.decouple_cell(nmt.compute_coupled_cell(f2,f2),cl_bias=clb22)
        np.savetxt(prefix+"_cl_%04d.txt"%i,
                   np.transpose([b.get_effective_ells(),
                                 cl22[0],cl22[1],cl22[2],cl22[3]]))
    cld=np.loadtxt(prefix+"_cl_%04d.txt"%i,unpack=True)
    cl22_all.append([cld[1],cld[2],cld[3],cld[4]])
cl22_all=np.array(cl22_all)

#Plot results
if plotres :
    cols=plt.cm.rainbow(np.linspace(0,1,6))
    plt.figure()
    plt.errorbar(b.get_effective_ells(),
                 np.mean(cl22_all,axis=0)[0]/cl22_th[0]-1,
                 yerr=np.std(cl22_all,axis=0)[0]/cl22_th[0]/np.sqrt(nsims+0.),
                 label='$EE$',fmt='bo')
    plt.errorbar(b.get_effective_ells(),
                 np.mean(cl22_all,axis=0)[3]/cl22_th[3]-1,
                 yerr=np.std(cl22_all,axis=0)[3]/cl22_th[3]/np.sqrt(nsims+0.),
                 label='$BB$',fmt='go')
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$\\Delta C_\\ell/C_\\ell$',fontsize=16)
    plt.legend(loc='lower right',frameon=False,fontsize=16)
    plt.xscale('log')

    ic=0
    plt.figure()
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[0],
             label='$EE$',c=cols[ic]);
    plt.plot(b.get_effective_ells(),cl22_th[0],'--',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[1],
             label='$EB$',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[3],
             label='$BB$',c=cols[ic]); ic+=1
    plt.loglog()
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='lower left',frameon=False,fontsize=14,ncol=2)
    plt.show()

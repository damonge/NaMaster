import numpy as np
import matplotlib.pyplot as plt
import os
import healpy as hp
import pymaster as nmt
import sys

DTOR=np.pi/180

if len(sys.argv)!=7 :
    print "python check_sph.py nside w_mask w_cont nsims plotres aposize "
    exit(1)

nside  =  int(sys.argv[1])
w_mask =  int(sys.argv[2])
w_cont =  int(sys.argv[3])
nsims  =  int(sys.argv[4])
plotres=  int(sys.argv[5])
aposize=float(sys.argv[6])

predir="tests_sph"
os.system("mkdir -p "+predir)
prefix=predir+"/run_ns%d_mask%d_cont%d_apo%.2lf"%(nside,w_mask,w_cont,aposize)
fname_mask=prefix+"_mask"

#This just generates the theory power spectra
if not os.path.isfile('cls_flat.txt') :
    import pyccl as ccl
    z=np.linspace(0.2,1.2,256)
    pz=np.exp(-((z-0.7)/0.1)**2)
    bz=np.ones_like(z)
    plt.figure(); plt.plot(z,pz); plt.xlabel('$z$',fontsize=16); plt.ylabel('$p(z)$',fontsize=16)
    cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=0.67, A_s=2.1e-9, n_s=0.96)
    clust=ccl.ClTracerNumberCounts(cosmo,False,False,z=z,n=pz,bias=bz)
    lens=ccl.ClTracerLensing(cosmo,False,z=z,n=pz)
    ell=np.arange(20000)
    cltt=ccl.angular_cl(cosmo,clust,clust,ell)
    clte=ccl.angular_cl(cosmo,clust,lens,ell)
    clee=ccl.angular_cl(cosmo,lens,lens,ell)
    clbb=np.zeros_like(clee)

    np.savetxt("cls_flat.txt",np.transpose([ell,cltt,clee,clbb,clte]))
ell,cltt,clee,clbb,clte=np.loadtxt("cls_flat.txt",unpack=True)
ell=ell[:3*nside]; cltt=cltt[:3*nside]; clee=clee[:3*nside]; clbb=clbb[:3*nside]; clte=clte[:3*nside];
cltt[0]=0; clee[0]=0; clbb[0]=0; clte[0]=0; 
if plotres :
    plt.figure()
    plt.plot(ell,cltt,'r-',label='$\\delta_g-\\delta_g$')
    plt.plot(ell,clte,'g-',label='$\\delta_g-\\gamma_E$')
    plt.plot(ell,clee,'b-',label='$\\gamma_E-\\gamma_E$')
    plt.loglog()
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='lower left',frameon=False,fontsize=16,labelspacing=0.1)

if w_cont :
    tilt_fg=-2.0
    l0_fg=100.
    clttfg=1E-5*((ell+10.)/(l0_fg+10.))**tilt_fg
    cleefg=1E-5*((ell+30.)/(l0_fg+30.))**tilt_fg
    cltefg=0.9*np.sqrt(clttfg*cleefg)
    clbbfg=cleefg
    clttfg[0]=0; cleefg[0]=0; clbbfg[0]=0; cltefg[0]=0; 
    if plotres :
        plt.plot(ell,clttfg,'r--',label='${\\rm FG},\\,TT$')
        plt.plot(ell,cltefg,'g--',label='${\\rm FG},\\,TE$')
        plt.plot(ell,cleefg,'b--',label='${\\rm FG},\\,EE$')
        plt.plot(ell,clbbfg,'y--',label='${\\rm FG},\\,BB$')

#This generates the mask with some padding and some holes
np.random.seed(1001)
fsky=0.1
rholes=1.
if not os.path.isfile(fname_mask+'.fits') :
    print "Generating mask"
    if w_mask :
        theta,phi=hp.pix2ang(nside,np.arange(hp.nside2npix(nside))); phi[phi>=np.pi]-=2*np.pi
        cth0=-np.sqrt(fsky); cthf= np.sqrt(fsky)
        theta0=np.arccos(cthf); thetaf=np.arccos(cth0)
        phi0=-np.pi*np.sqrt(fsky); phif=np.pi*np.sqrt(fsky)
        ids=np.where((theta>theta0) & (theta<thetaf) &
                     (phi>phi0) & (phi<phif))[0]
        mask_raw=np.zeros(hp.nside2npix(nside)); mask_raw[ids]=1.
        nholes=15
        cths=cth0+(cthf-cth0)*np.random.rand(nholes)
        phis=phi0+(phif-phi0)*np.random.rand(nholes)
        ths=np.arccos(cths)
        vs=np.transpose(np.array([np.sin(ths)*np.cos(phis),np.sin(ths)*np.sin(phis),np.cos(ths)]))
        for i in np.arange(nholes) :
            v=vs[i]
            mask_raw[hp.query_disc(nside,vs[i],rholes*np.pi/180)]=0
        if aposize>0 :
            mask=nmt.mask_apodization(mask_raw,aposize,apotype='C1')
        else :
            mask=mask_raw
    else :
        mask=np.ones(hp.nside2npix(nside))

    hp.write_map(fname_mask+".fits",mask)
mask=hp.read_map(fname_mask+".fits")
if plotres :
    hp.mollview(mask)

if w_cont :
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
alpha_cont_0=0.1
alpha_cont_2=0.1
def get_fields() :
    mppt,mppq,mppu=nmt.synfast_spherical(nside,[cltt,clee,clbb,clte],pol=True)
    if w_cont :
        mppt+=alpha_cont_0*fgt
        mppq+=alpha_cont_2*fgq
        mppu+=alpha_cont_2*fgu
        ff0=nmt.NmtField(mask,[mppt],templates=[[fgt]])
        ff2=nmt.NmtField(mask,[mppq,mppu],[[fgq,fgu]])
    else :
        ff0=nmt.NmtField(mask,[mppt])
        ff2=nmt.NmtField(mask,[mppq,mppu])
    return mppt,mppq,mppu,ff0,ff2
mpt,mpq,mpu,f0,f2=get_fields()
    
if plotres :
    hp.mollview((mpt*mask).flatten(),title='$\\delta_g$')
    hp.mollview((mpq*mask).flatten(),title='$\\gamma_1$')
    hp.mollview((mpu*mask).flatten(),title='$\\gamma_2$')

#Compute deprojection bias
if w_cont :
    clb00=nmt.deprojection_bias(f0,f0,[cltt])
    clb02=nmt.deprojection_bias(f0,f2,[clte,0*clte])
    clb22=nmt.deprojection_bias(f2,f2,[clee,0*clee,0*clbb,clbb])
else :
    clb00=None;
    clb02=None;
    clb22=None;
    
#Use initial fields to generate coupling matrix
w00=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w00.dat") :
    print "Computing 00"
    w00.compute_coupling_matrix(f0,f0,b)
    w00.write_to(prefix+"_w00.dat");
else :
    w00.read_from(prefix+"_w00.dat")
w02=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w02.dat") :
    print "Computing 02"
    w02.compute_coupling_matrix(f0,f2,b)
    w02.write_to(prefix+"_w02.dat");
else :
    w02.read_from(prefix+"_w02.dat")
w22=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w22.dat") :
    print "Computing 22"
    w22.compute_coupling_matrix(f2,f2,b)
    w22.write_to(prefix+"_w22.dat");
else :
    w22.read_from(prefix+"_w22.dat")

#Generate theory prediction
cl00_th=w00.decouple_cell(w00.couple_cell([cltt]))
cl02_th=w02.decouple_cell(w02.couple_cell([clte,0*clte]))
cl22_th=w22.decouple_cell(w22.couple_cell([clee,0*clee,0*clbb,clbb]))
np.savetxt(prefix+"_cl_th.txt",
           np.transpose([b.get_effective_ells(),cl00_th[0],cl02_th[0],cl02_th[1],
                         cl22_th[0],cl22_th[1],cl22_th[2],cl22_th[3]]))

#Compute mean and variance over nsims simulations
cl00_all=[]
cl02_all=[]
cl22_all=[]
for i in np.arange(nsims) :
    if i%100==0 :
        print "%d-th sim"%i

    if not os.path.isfile(prefix+"_cl_%04d.txt"%i) :
        mpt,mpq,mpu,f0,f2=get_fields()
        cl00=w00.decouple_cell(nmt.compute_coupled_cell(f0,f0),cl_bias=clb00)
        cl02=w02.decouple_cell(nmt.compute_coupled_cell(f0,f2),cl_bias=clb02)
        cl22=w22.decouple_cell(nmt.compute_coupled_cell(f2,f2),cl_bias=clb22)
        np.savetxt(prefix+"_cl_%04d.txt"%i,
                   np.transpose([b.get_effective_ells(),cl00[0],cl02[0],cl02[1],
                                 cl22[0],cl22[1],cl22[2],cl22[3]]))
    cld=np.loadtxt(prefix+"_cl_%04d.txt"%i,unpack=True)
    cl00_all.append([cld[1]])
    cl02_all.append([cld[2],cld[3]])
    cl22_all.append([cld[4],cld[5],cld[6],cld[7]])
cl00_all=np.array(cl00_all)
cl02_all=np.array(cl02_all)
cl22_all=np.array(cl22_all)

#Plot results
if plotres :
    cols=plt.cm.rainbow(np.linspace(0,1,6))
    plt.figure()
    plt.errorbar(b.get_effective_ells(),np.mean(cl00_all,axis=0)[0]/cl00_th[0]-1,yerr=np.std(cl00_all,axis=0)[0]/cl00_th[0]/np.sqrt(nsims+0.),label='$\\delta_g-\\delta_g$',fmt='ro')
    plt.errorbar(b.get_effective_ells(),np.mean(cl02_all,axis=0)[0]/cl02_th[0]-1,yerr=np.std(cl02_all,axis=0)[0]/cl02_th[0]/np.sqrt(nsims+0.),label='$\\delta_g-\\gamma_E$',fmt='go')
    plt.errorbar(b.get_effective_ells(),np.mean(cl22_all,axis=0)[0]/cl22_th[0]-1,yerr=np.std(cl22_all,axis=0)[0]/cl22_th[0]/np.sqrt(nsims+0.),label='$\\gamma_E-\\gamma_E$',fmt='bo')
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$\\Delta C_\\ell/C_\\ell$',fontsize=16)
    plt.legend(loc='lower right',frameon=False,fontsize=16)
    plt.xscale('log')

    ic=0
    plt.figure()
    plt.plot(b.get_effective_ells(),np.mean(cl00_all,axis=0)[0],
             label='$\\delta_g-\\delta_g$',c=cols[ic])
    plt.plot(b.get_effective_ells(),cl00_th[0],'--',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl02_all,axis=0)[0],
             label='$\\delta_g-\\gamma_E$',c=cols[ic]);
    plt.plot(b.get_effective_ells(),cl02_th[0],'--',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl02_all,axis=0)[1],
             label='$\\delta_g-\\gamma_B$',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[0],
             label='$\\gamma_E-\\gamma_E$',c=cols[ic]);
    plt.plot(b.get_effective_ells(),cl22_th[0],'--',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[1],
             label='$\\gamma_E-\\gamma_B$',c=cols[ic]); ic+=1
    plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[3],
             label='$\\gamma_B-\\gamma_B$',c=cols[ic]); ic+=1
    plt.loglog()
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='lower left',frameon=False,fontsize=14,ncol=2)
    plt.show()

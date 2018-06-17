import numpy as np
import matplotlib.pyplot as plt
import os
import healpy as hp
import pymaster as nmt
import sys
import time
DTOR=np.pi/180

if len(sys.argv)!=7 :
    print "python check_sph.py nside w_mask w_cont nsims plotres aposize "
    exit(1)
t0 = time.time()
nside  =  int(sys.argv[1])
w_mask =  int(sys.argv[2])
w_cont =  int(sys.argv[3])
nsims  =  int(sys.argv[4])
plotres=  int(sys.argv[5])
aposize=float(sys.argv[6])
star_map = np.load('/global/cscratch1/sd/jsanch87/tests_sph/catsim_nside512_stars.npy')
star_mask = np.load('/global/cscratch1/sd/jsanch87/tests_sph/catsim_nside512_mask.npy')
footprint = np.load('/project/projectdirs/lsst/old/LSSWG/DepthVariations/coaddM5Data_masked_rBand_RepulsiveRandomDitherFieldPerVisit.npz')
mask = np.logical_not(footprint['mask'])
mask = hp.ud_grade(mask,nside_out=nside)
star_map = hp.ud_grade(star_map,nside_out=nside)
delta_star = np.zeros_like(star_map)
delta_star[mask] = star_map[mask]/np.mean(star_map[mask])-1
alpha_star = 0.03 # We fix the stellar contamination at 3 percent level
K = alpha_star/np.sqrt(np.mean(delta_star**2))
print 'K=',K
cl_stars = hp.anafast(K*delta_star)
ells_stars = np.arange(len(cl_stars))
shot_noise_wgt = np.zeros(12*nside**2) # We are going to add shot noise proportional to the stellar density
# This is because the larger number, the smaller the number of detected galaxies will be.
shot_noise_wgt[mask] = star_map[mask]/np.max(star_map[mask])
n_bar = 1./(1./3600.*np.pi/180*np.pi/180)*46 # Expected number density for LSST (we assume ~46 galaxies/sq-arcmin and convert to  galaxies/sterad)
base_sn = 1./n_bar
print 'Shot noise level', base_sn
predir="/global/cscratch1/sd/jsanch87/tests_NaMaster_clustering"
#os.system("mkdir -p "+predir)
prefix=predir+"/run_ns%d_mask%d_cont%d_apo%.2lf"%(nside,w_mask,w_cont,aposize)
fname_mask=prefix+"_mask"

#This just generates the theory power spectra
if not os.path.isfile('cls_clustering.txt') :
    import pyccl as ccl
    z=np.linspace(0.2,1.2,256)
    pz=np.exp(-0.5*((z-0.5)/0.075)**2)
    bz=np.ones_like(z)
    plt.figure(); plt.plot(z,pz); plt.xlabel('$z$',fontsize=16); plt.ylabel('$p(z)$',fontsize=16)
    cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=0.67, A_s=2.1e-9, n_s=0.96)
    clust=ccl.ClTracerNumberCounts(cosmo,False,False,z=z,n=pz,bias=bz)
    lens=ccl.ClTracerLensing(cosmo,False,z=z,n=pz)
    ell=np.arange(40000)
    cltt=ccl.angular_cl(cosmo,clust,clust,ell)
    clte=ccl.angular_cl(cosmo,clust,lens,ell)
    clee=ccl.angular_cl(cosmo,lens,lens,ell)
    clbb=np.zeros_like(clee)

    np.savetxt("cls_clustering.txt",np.transpose([ell,cltt,clee,clbb,clte]))
ell,cltt,clee,clbb,clte=np.loadtxt("cls_clustering.txt",unpack=True)
ell=ell[:3*nside]; cltt=cltt[:3*nside]; clee=clee[:3*nside]; clbb=clbb[:3*nside]; clte=clte[:3*nside];
cltt[0]=0; clee[0]=0; clbb[0]=0; clte[0]=0; 
if plotres :
    plt.figure()
    plt.plot(ell,cltt,'r-',label='$\\delta_g-\\delta_g$')
    plt.plot(ell,clte,'g-',label='$\\delta_g-\\gamma_E$')
    plt.plot(ell,clee,'b-',label='$\\gamma_E-\\gamma_E$')
    if w_cont:
       plt.plot(ells_stars,cl_stars, label='$K(\\delta_*-\\delta_*)$')
    plt.loglog()
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='lower left',frameon=False,fontsize=16,labelspacing=0.1)
    plt.savefig(predir+"/Cls_input_ns%d_mask%d_cont%d_apo%.2lf.png"%(nside,w_mask,w_cont,aposize))


#This generates the mask with some padding and some holes
np.random.seed(1001)
rholes=0.05
if not os.path.isfile(fname_mask+'.fits') :
    print "Generating mask"
    if w_mask :  
        mask_raw=mask
        nholes=15
        cths=-1+2*np.random.rand(nholes)
        phis=0+2*np.pi*np.random.rand(nholes)
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
fsky = np.mean(mask)
if plotres :
    hp.mollview(mask)
print('fsky: ', np.mean(mask))
if w_cont :
    if not os.path.isfile(prefix+"_contaminants.fits") :
        #fgt,fgq,fgu=nmt.synfast_spherical(nside,[clttfg,cleefg,clbbfg,cltefg],pol=True)
        shot_noise_map = nmt.synfast_spherical(nside,base_sn*np.ones(3*nside))[0]
        plt.figure()
        plt.plot(hp.anafast(shot_noise_map))
        #shot_noise_map = shot_noise_map-np.min(shot_noise_map)
        #hp.mollview(shot_noise_map)
        plt.show()
        hp.mollview(shot_noise_map*shot_noise_wgt)
        plt.show()
        hp.mollview(K*delta_star)
        plt.show()
        fgt =  K*delta_star + shot_noise_map*shot_noise_wgt
        fgq = shot_noise_map*shot_noise_wgt
        fgu = fgq
        hp.write_map(prefix+"_contaminants.fits",[fgt,fgq,fgu])
    else :
        fgt,fgq,fgu=hp.read_map(prefix+"_contaminants.fits",field=[0,1,2],verbose=False)

#Binning scheme
d_ell=int(1./fsky)
t0 = time.time()
b=nmt.NmtBin(nside,nlb=d_ell)
t1 = time.time()
print('Binning time: ', t1-t0)
#Generate some initial fields
print " - Res: %.3lf arcmin. "%(np.sqrt(4*np.pi*(180*60/np.pi)**2/hp.nside2npix(nside)))
def get_fields() :
    mppt,mppq,mppu=nmt.synfast_spherical(nside,[cltt,clee,clbb,clte],pol=True)
    if w_cont : 
        mppt+=fgt
        mppq+=fgq
        mppu+=fgu
        ff0=nmt.NmtField(mask,[mppt],templates=[[fgt]])
        ff2=nmt.NmtField(mask,[mppq,mppu],[[fgq,fgu]])
        #ff0=nmt.NmtField(mask,[mppt])
        #ff2=nmt.NmtField(mask,[mppq,mppu])
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
t0 = time.time()    
#Use initial fields to generate coupling matrix
w00=nmt.NmtWorkspace();
t1 = time.time()
print('Workspace generation time', t1-t0)
if not os.path.isfile(prefix+"_w00.dat") :
    print "Computing 00"
    t0 = time.time()
    w00.compute_coupling_matrix(f0,f0,b)
    t1 = time.time()
    print('Coupling matrix (w00) time', t1-t0)
    w00.write_to(prefix+"_w00.dat");
else :
    w00.read_from(prefix+"_w00.dat")
w02=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w02.dat") :
    print "Computing 02"
    t0 = time.time()
    w02.compute_coupling_matrix(f0,f2,b)
    t1 = time.time()
    print('Coupling matrix (w02) time', t1-t0)
    w02.write_to(prefix+"_w02.dat");
else :
    w02.read_from(prefix+"_w02.dat")
w22=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w22.dat") :
    print "Computing 22"
    t0 = time.time()
    w22.compute_coupling_matrix(f2,f2,b)
    t1 = time.time()
    print('Coupling matrix (w22) time', t1-t0)
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
        t0 = time.time()
        cl00=w00.decouple_cell(nmt.compute_coupled_cell(f0,f0),cl_bias=clb00)
        t1 = time.time()
        print('Time to decouple cl00', t1-t0)
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
    plt.savefig(predir+"/Cls_ratio_ns%d_mask%d_cont%d_apo%.2lf.png"%(nside,w_mask,w_cont,aposize))
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
    plt.savefig(predir+"/Cls_output_ns%d_mask%d_cont%d_apo%.2lf.png"%(nside,w_mask,w_cont,aposize))
    plt.show()

#cw=nmt.NmtCovarianceWorkspace()
#cw.compute_coupling_coefficients(w00,w00) #<- This is the time-consuming operation
#covar=nmt.gaussian_covariance(cw,cltt,cltt,cltt,cltt)
#err_th = cl00_th[0]/np.sqrt(d_ell*fsky)*np.sqrt(2./(2*b.get_effective_ells()+1.))
#plt.figure()
#plt.plot(b.get_effective_ells(),np.std(cl00_all,axis=0)[0]/(np.sqrt(np.diag(covar))),label='Gaussian')
#plt.plot(b.get_effective_ells(),np.std(cl00_all,axis=0)[0]/err_th,label='Theory')
#plt.grid()
#plt.xlabel(r'$\ell$')
#plt.ylabel(r'$\Delta C_{\ell,disp}/\Delta C_{\ell,th}$')
#plt.xscale('log')
#plt.legend(loc='best')
#plt.show()

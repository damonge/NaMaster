import numpy as np
import matplotlib.pyplot as plt
import os
import healpy as hp
import pymaster as nmt
import sys

DTOR=np.pi/180

if len(sys.argv)!=8 :
    print "python check_sph.py nside w_mask w_cont w_nvar nsims plotres aposize "
    exit(1)

nside  =  int(sys.argv[1])
w_mask =  int(sys.argv[2])
w_cont =  int(sys.argv[3])
w_nvar =  int(sys.argv[4])
nsims  =  int(sys.argv[5])
plotres=  int(sys.argv[6])
aposize=float(sys.argv[7])

#Switch off contaminants if there's no mask
if not w_mask :
    w_cont=False
    w_nvar=False
    
predir="tests_sph"
os.system("mkdir -p "+predir)
prefix=predir+"/run_ns%d_mask%d_cont%d_apo%.2lf"%(nside,w_mask,w_cont,aposize)
fname_mask=prefix+"_mask"

l,cltt,clee,clbb,clte,nltt,nlee,nlbb,nlte=np.loadtxt("data/cls_lss.txt",unpack=True)
cltt=cltt[:3*nside]; clee=clee[:3*nside]; clbb=clbb[:3*nside]; clte=clte[:3*nside]; 
nltt=nltt[:3*nside]; nlee=nlee[:3*nside]; nlbb=nlbb[:3*nside]; nlte=nlte[:3*nside]; 
np.random.seed(1001)


if not os.path.isfile(fname_mask+'.fits') :
    if w_nvar :
        depth_nvar=hp.read_map("data/cont_lss_nvar_ns%d.fits"%nside,verbose=False)
        depth_nvar[depth_nvar<0.8]=0
    else :
        depth_nvar=np.ones(hp.nside2npix(nside))
    
    if w_mask :
        depth_ivar=np.zeros_like(depth_nvar); depth_ivar[depth_nvar>0.1]=1./depth_nvar[depth_nvar>0.1]
        mask_raw=hp.read_map("data/mask_lss_ns%d.fits"%nside,verbose=False)
        if aposize>0 :
            mask=nmt.mask_apodization(mask_raw,aposize,apotype='C1')
        else :
            mask=mask_raw
        mask*=depth_ivar
    else :
        mask=np.ones(hp.nside2npix(nside))
    
    depth_nvar=np.sqrt(depth_nvar)
    hp.write_map(fname_mask+".fits",mask,overwrite=True)
mask=hp.read_map(fname_mask+".fits",verbose=False)
if plotres :
    hp.mollview(mask)
fsky=np.mean(mask/np.amax(mask));

if w_cont :
    fgt=np.zeros_like(mask); fgq=np.zeros_like(mask); fgu=np.zeros_like(mask)
    t=hp.read_map("data/cont_lss_star_ns%d.fits"%nside,verbose=False); fgt+=t
    t=hp.read_map("data/cont_lss_dust_ns%d.fits"%nside,verbose=False); fgt+=t
    q,u=hp.read_map("data/cont_wl_psf_ns%d.fits"%nside,field=[0,1],verbose=False); fgq+=q; fgu+=u;
    q,u=hp.read_map("data/cont_wl_ss_ns%d.fits"%nside,field=[0,1],verbose=False); fgq+=q; fgu+=u;

#Binning scheme
d_ell=int(1./fsky)
b=nmt.NmtBin(nside,nlb=d_ell)

if plotres :
    hp.mollview(fgt*mask)
    hp.mollview(fgq*mask)
    hp.mollview(fgu*mask)

#Generate some initial fields
print " - Res: %.3lf arcmin. "%(np.sqrt(4*np.pi*(180*60/np.pi)**2/hp.nside2npix(nside)))
def get_fields() :
    st,sq,su=nmt.synfast_spherical(nside,[cltt,clee,clbb,clte],pol=True)
    nt,nq,nu=nmt.synfast_spherical(nside,[nltt,nlee,nlbb,nlte],pol=True)
    st+=depth_nvar*nt
    sq+=depth_nvar*nq
    su+=depth_nvar*nu
    if w_cont :
        st+=fgt
        sq+=fgq
        su+=fgu
        ff0=nmt.NmtField(mask,[st],templates=[[fgt]])
        ff2=nmt.NmtField(mask,[sq,su],[[fgq,fgu]])
    else :
        ff0=nmt.NmtField(mask,[st])
        ff2=nmt.NmtField(mask,[sq,su])
    return st,sq,su,ff0,ff2
mpt,mpq,mpu,f0,f2=get_fields()
    
if plotres :
    hp.mollview((mpt*mask).flatten(),title='$\\delta_g$')
    hp.mollview((mpq*mask).flatten(),title='$\\gamma_1$')
    hp.mollview((mpu*mask).flatten(),title='$\\gamma_2$')
plt.show(); exit(1)

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

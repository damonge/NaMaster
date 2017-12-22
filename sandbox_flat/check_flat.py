import numpy as np
import matplotlib.pyplot as plt
import os
import flatmaps as fm
import pymaster as nmt

DTOR=np.pi/180
fname_mask="mask_flat"

#Switch both to false if you want a simple mask
make_holes=True
make_padding=True
if not make_holes :
    fname_mask+="_no_holes"
if not make_padding :
    fname_mask+="_no_padding"

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
plt.figure()
plt.plot(ell,cltt,'r-',label='$\\delta_g-\\delta_g$')
plt.plot(ell,clte,'g-',label='$\\delta_g-\\gamma_E$')
plt.plot(ell,clee,'b-',label='$\\gamma_E-\\gamma_E$')
plt.loglog()
plt.xlabel('$\\ell$',fontsize=16)
plt.ylabel('$C_\ell$',fontsize=16)
plt.legend(loc='lower left',frameon=False,fontsize=16,labelspacing=0.1)

#This generates the mask with some padding and some holes
np.random.seed(1001)
if not os.path.isfile(fname_mask+'.npz') :
    mi=fm.FlatMapInfo([212.5,222.],[-2.,2.],nx=792,ny=336)
    print "Generating mask"
    mask=np.ones([mi.ny,mi.nx])

    if make_padding :
        mask[:mi.ny/7 ,:]=0; mask[ 6*mi.ny/7: ,:]=0
        mask[:,:mi.nx/16]=0; mask[:,15*mi.nx/16:]=0
        
    mask=mask.flatten()

    if make_holes :
        xarr=np.ones(mi.ny)[:,None]*np.arange(mi.nx)[None,:]*mi.lx/mi.nx
        yarr=np.ones(mi.nx)[None,:]*np.arange(mi.ny)[:,None]*mi.ly/mi.ny
        def dig_hole(x,y,r) :
            rad=(np.sqrt((xarr-x)**2+(yarr-y)**2)).flatten()
            return np.where(rad<r)[0]
        nholes=50
        for i in np.arange(nholes) :
            r=np.random.rand(3)
            mask[dig_hole(r[0]*mi.lx,r[1]*mi.ly,(0.005+0.02*r[2])*np.sqrt(mi.lx*mi.ly))]=0.
            
    mi.write_flat_map(fname_mask,mask)
mi,mask=fm.read_flat_map(fname_mask+'.npz')
mi=fm.FlatMapInfo([212.5,222.],[-2.,2.],nx=792,ny=336)
mask=mask.reshape([mi.ny,mi.nx])
mi.view_map(mask.flatten(),addColorbar=False,colorMax=1,colorMin=0)

#Binning scheme
lx_rad=mi.lx*DTOR
ly_rad=mi.ly*DTOR
ell_min=max(2*np.pi/(mi.lx*DTOR),2*np.pi/(mi.ly*DTOR))
ell_max=min(mi.nx*np.pi/(mi.lx*DTOR),mi.ny*np.pi/(mi.ly*DTOR))
d_ell=2*ell_min
n_ell=int((ell_max-ell_min)/d_ell)-1
l_bpw=np.zeros([2,n_ell])
l_bpw[0,:]=ell_min+np.arange(n_ell)*d_ell
l_bpw[1,:]=l_bpw[0,:]+d_ell
b=nmt.NmtBinFlat(l_bpw[0,:],l_bpw[1,:])

#Generate some initial fields
print " - Res(x): %.3lf arcmin. Res(y): %.3lf arcmin."%(mi.lx*60/mi.nx,mi.ly*60/mi.ny)
print " - lmax = %d, lmin = %d"%(int(ell_max),int(ell_min))
mpt,mpq,mpu=nmt.synfast_flat(int(mi.nx),int(mi.ny),mi.lx*DTOR,mi.ly*DTOR,[cltt,clee,clbb,clte],pol=True)
f0=nmt.NmtFieldFlat(mi.lx*DTOR,mi.ly*DTOR,mask,[mpt])
f2=nmt.NmtFieldFlat(mi.lx*DTOR,mi.ly*DTOR,mask,[mpq,mpu])
mi.view_map((mpt*mask).flatten(),title='$\\delta_g$',addColorbar=False)
mi.view_map((mpq*mask).flatten(),title='$\\gamma_1$',addColorbar=False)
mi.view_map((mpu*mask).flatten(),title='$\\gamma_2$',addColorbar=False)

#Use initial fields to generate coupling matrix
w00=nmt.NmtWorkspaceFlat();
if not os.path.isfile("w00_"+fname_mask+".dat") :
    print "Computing 00"
    w00.compute_coupling_matrix(f0,f0,b)
    w00.write_to("w00_"+fname_mask+".dat");
else :
    w00.read_from("w00_"+fname_mask+".dat")
w02=nmt.NmtWorkspaceFlat();
if not os.path.isfile("w02_"+fname_mask+".dat") :
    print "Computing 02"
    w02.compute_coupling_matrix(f0,f2,b)
    w02.write_to("w02_"+fname_mask+".dat");
else :
    w02.read_from("w02_"+fname_mask+".dat")
w22=nmt.NmtWorkspaceFlat();
if not os.path.isfile("w22_"+fname_mask+".dat") :
    print "Computing 22"
    w22.compute_coupling_matrix(f2,f2,b)
    w22.write_to("w22_"+fname_mask+".dat");
else :
    w22.read_from("w22_"+fname_mask+".dat")

#Generate theory prediction
cl00_th=w00.decouple_cell(w00.couple_cell(ell,np.array([cltt])))
cl02_th=w02.decouple_cell(w02.couple_cell(ell,np.array([clte,0*clte])))
cl22_th=w22.decouple_cell(w22.couple_cell(ell,np.array([clee,0*clee,0*clbb,clbb])))

#Compute mean and variance over nsims simulations
nsims=1000
cl00_all=[]
cl02_all=[]
cl22_all=[]
for i in np.arange(nsims) :
    mpt,mpq,mpu=nmt.synfast_flat(mi.nx,mi.ny,mi.lx*DTOR,mi.ly*DTOR,[cltt,clee,clbb,clte],pol=True)
    if i%100==0 :
        print "%d-th sim"%i
    f0=nmt.NmtFieldFlat(mi.lx*DTOR,mi.ly*DTOR,mask,[mpt])
    f2=nmt.NmtFieldFlat(mi.lx*DTOR,mi.ly*DTOR,mask,[mpq,mpu])
    cl00_all.append(w00.decouple_cell(nmt.compute_coupled_cell_flat(f0,f0,b)))
    cl02_all.append(w02.decouple_cell(nmt.compute_coupled_cell_flat(f0,f2,b)))
    cl22_all.append(w22.decouple_cell(nmt.compute_coupled_cell_flat(f2,f2,b)))
cl00_all=np.array(cl00_all)
cl02_all=np.array(cl02_all)
cl22_all=np.array(cl22_all)

#Plot results
cols=plt.cm.rainbow(np.linspace(0,1,6))
ic=0
plt.figure()
plt.errorbar(b.get_effective_ells(),np.mean(cl00_all,axis=0)[0]/cl00_th[0]-1,yerr=np.std(cl00_all,axis=0)[0]/cl00_th[0]/np.sqrt(nsims+0.),label='$\\delta_g-\\delta_g$',fmt='ro')
plt.errorbar(b.get_effective_ells(),np.mean(cl02_all,axis=0)[0]/cl02_th[0]-1,yerr=np.std(cl02_all,axis=0)[0]/cl02_th[0]/np.sqrt(nsims+0.),label='$\\delta_g-\\gamma_E$',fmt='go')
plt.errorbar(b.get_effective_ells(),np.mean(cl22_all,axis=0)[0]/cl22_th[0]-1,yerr=np.std(cl22_all,axis=0)[0]/cl22_th[0]/np.sqrt(nsims+0.),label='$\\gamma_E-\\gamma_E$',fmt='bo')
plt.xlabel('$\\ell$',fontsize=16)
plt.ylabel('$C_\ell$',fontsize=16)
plt.legend(loc='lower right',frameon=False,fontsize=16)
plt.xscale('log')

plt.figure()
plt.plot(b.get_effective_ells(),np.mean(cl00_all,axis=0)[0],label='$\\delta_g-\\delta_g$',c=cols[ic])
plt.plot(b.get_effective_ells(),cl00_th[0],'--',c=cols[ic]); ic+=1
plt.plot(b.get_effective_ells(),np.mean(cl02_all,axis=0)[0],label='$\\delta_g-\\gamma_E$',c=cols[ic]);
plt.plot(b.get_effective_ells(),cl02_th[0],'--',c=cols[ic]); ic+=1
plt.plot(b.get_effective_ells(),np.mean(cl02_all,axis=0)[1],label='$\\delta_g-\\gamma_B$',c=cols[ic]); ic+=1
plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[0],label='$\\gamma_E-\\gamma_E$',c=cols[ic]);
plt.plot(b.get_effective_ells(),cl22_th[0],'--',c=cols[ic]); ic+=1
plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[1],label='$\\gamma_E-\\gamma_B$',c=cols[ic]); ic+=1
plt.plot(b.get_effective_ells(),np.mean(cl22_all,axis=0)[3],label='$\\gamma_B-\\gamma_B$',c=cols[ic]); ic+=1
plt.loglog()
plt.xlabel('$\\ell$',fontsize=16)
plt.ylabel('$C_\ell$',fontsize=16)
plt.legend(loc='lower left',frameon=False,fontsize=14,ncol=2)
plt.show()

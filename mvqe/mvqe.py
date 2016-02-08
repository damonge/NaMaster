import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cg_invert as inv
import scipy.linalg as llg

nside=64#
lmax=3*nside-1
npix=hp.nside2npix(nside)
sigma_n=0.3
glb_niter=0

def sl_fun(ell) :
    return (ell/60.)**(-1.0)

def nl_fun(ell) :
    return sigma_n*np.ones_like(ell)

lcut=100

l_arr=np.arange(lmax+1)
sl_arr=sl_fun(l_arr); sl_arr[0]=sl_arr[1]
nl_arr=nl_fun(l_arr)
sl_prop_arr=np.array(sl_arr); sl_prop_arr[l_arr>lcut]=sl_arr[lcut]
isl_prop_arr=1./sl_prop_arr

mask=np.ones(npix)
mask=hp.ud_grade(hp.read_map("mask_2mass.fits",verbose=False),nside)
ipix_seen=np.where(mask>0.1)[0]
npix_seen=len(ipix_seen)
map_s=hp.synfast(sl_arr,nside,verbose=False)
map_n=hp.synfast(nl_arr,nside,verbose=False)
sigma_noise=np.sqrt(sigma_n*npix/(4*np.pi))
map_d=(map_s+map_n)*mask
#hp.mollview(map_d); plt.show()

plt.plot(l_arr,sl_arr)
plt.plot(l_arr,nl_arr)
plt.plot(l_arr,sl_prop_arr)
plt.plot(l_arr,sl_arr+nl_arr)
#plt.plot(l_arr,hp.anafast(map_n))
#plt.plot(l_arr,hp.anafast(map_s))
#plt.plot(l_arr,hp.anafast(map_d))
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.show()

res_map=hp.alm2map(hp.map2alm(map_d,iter=0),nside,verbose=False)-map_d
res_limit=np.sqrt(np.sum(res_map**2)/np.sum((map_s+map_n)**2))
print res_limit
hp.mollview(res_map); plt.show()

nbl=4
nb=(lmax-2)/nbl
bins=2+np.arange(nb+1)*nbl
print bins

pars=inv.ParamMVQE(ipix_seen,sl_arr,sigma_noise,nside,lcut,1E-3,True,bins)
corr_cb=inv.correlated_cl(pars,map_d)
pars.plot_stuff=False
fisher=inv.fisher_avg(pars,100)
uncorr_cb=llg.solve(fisher,corr_cb)
lbin=0.5*(bins[:-1]+bins[1:])
np.savetxt("uncorrs.txt",np.transpose([lbin,uncorr_cb,corr_cb]))
lbin,uncorr_cb,corr_cb=np.loadtxt("uncorrs.txt",unpack=True)
plt.plot(lbin,uncorr_cb)
plt.plot(l_arr,sl_arr+nl_arr)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.show()

#inv.invert_covar(map_d,ipix_seen,sl_arr,sigma_noise,nside,lcut,1E-3,True)
#inv.invert_covar(map_d,ipix_seen,sl_arr,sigma_noise,nside,isl_prop_arr,res_limit*0.001,True)

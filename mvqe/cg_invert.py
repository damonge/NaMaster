import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg,LinearOperator

class ParamMVQE :
    ip_seen=None
    signal_cl=None
    sigma_noise=None
    nside=None
    lmax_precon=None
    cgtol=None
    plot_stuff=None
    bins=None
    cg_niter=0

    def __init__(self,ip_seen,signal_cl,sigma_noise,nside,lmax_precon,cgtol,plot_stuff,bins) :
        self.ip_seen=ip_seen
        self.signal_cl=signal_cl
        self.sigma_noise=sigma_noise
        self.nside=nside
        self.lmax_precon=lmax_precon
        self.cgtol=cgtol
        self.plot_stuff=plot_stuff
        self.bins=bins


def invert_covar(par,mp) :
    par.cg_niter=0
    npx_seen=len(par.ip_seen)
    npx=hp.nside2npix(par.nside)
    pixsize=4*np.pi/npx
    glb_niter=0
    
    cl_signal_precon=pixsize**2*(1./par.signal_cl[:par.lmax_precon]-1./par.signal_cl[par.lmax_precon])
    sig_precon=np.sqrt(pixsize/par.signal_cl[par.lmax_precon])

    def linop_covar(m) :
        m_full=np.zeros(npx); m_full[par.ip_seen]=m
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0),par.signal_cl/pixsize),
                         par.nside,verbose=False)+par.sigma_noise**2*m_full
        return m_out[par.ip_seen]

    def linop_precon(m) :
        m_full=np.zeros(npx); m_full[par.ip_seen]=m
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0,lmax=par.lmax_precon-1),
                                   cl_signal_precon),par.nside,verbose=False)/pixsize+sig_precon**2*m_full
        return m_out[par.ip_seen]

    def callback(m) :
        par.cg_niter+=1

    A=LinearOperator((npx_seen,npx_seen),matvec=linop_covar)
    M=LinearOperator((npx_seen,npx_seen),matvec=linop_precon)
    b=mp[par.ip_seen]

    m_solve,info=cg(A,b,x0=linop_precon(b),M=M,tol=par.cgtol,callback=callback)
    mp_out=np.zeros_like(mp); mp_out[par.ip_seen]=m_solve

    if par.plot_stuff :
        print "  Inverted in %d iterations"%(par.cg_niter)
        mp_fill=np.zeros_like(mp); mp_fill[par.ip_seen]=linop_covar(m_solve)
        hp.mollview(mp_fill)
        hp.mollview(mp)
        hp.mollview((mp_fill-mp)/np.std(mp));
        hp.mollview(mp_out);
        plt.show()

    return mp_out

def random_unit_variance(n_elements,do_gauss) :
    if do_gauss :
        return np.random.randn(n_elements)
    else :
        vr=np.random.rand(n_elements)
        v_out=np.ones(n_elements)
        v_out[np.where(vr<0.5)]=-1
        return v_out

def fisher_single(par,v) :
    # v -> seen pixels
    nb=len(par.bins)-1
    npix=hp.nside2npix(par.nside)
    npix_seen=len(par.ip_seen)
    lmax=3*par.nside-1
    larr=np.arange(lmax+1)
    fisher=np.zeros([nb,nb])
    pixsize=4*np.pi/hp.nside2npix(par.nside)
    
    v_map=np.zeros(npix); v_map[par.ip_seen]=v
    vcm1=invert_covar(par,v_map)
    v_lm=hp.map2alm(v_map,iter=0)
    vcm1_lm=hp.map2alm(vcm1,iter=0)
    for iba in np.arange(nb) :
#        print " Row %d"%iba
        transfer=np.zeros(lmax+1); transfer[par.bins[iba]:par.bins[iba+1]]=1.
        v_map2=hp.alm2map(hp.almxfl(v_lm,transfer),par.nside,verbose=False)/pixsize #Q_a * v
        v_map2cm1=invert_covar(par,v_map2) #C^-1 * Q_a * v
        va_lm=hp.map2alm(v_map2cm1,iter=0)
        cl_vcm1_va=(2*larr+1)*hp.alm2cl(vcm1_lm,alms2=va_lm)
        for ibb in np.arange(nb-iba)+iba :
            fisher[iba,ibb]=np.sum(cl_vcm1_va[par.bins[ibb]:par.bins[ibb+1]])/pixsize**2
            if iba!=ibb :
                fisher[ibb,iba]=fisher[iba,ibb]

    return fisher

def fisher_avg(par,n_vec) :
    nb=len(par.bins)-1
    npix_seen=len(par.ip_seen)
    fisher=np.zeros([nb,nb])
    
    for iv in np.arange(n_vec) :
        print "Realization #%d"%iv
        v=random_unit_variance(npix_seen,False)
        fisher+=fisher_single(par,v)
        mean_here=np.sum(fisher/(iv+1))
        print " Mean: %lE"%mean_here
        
    fisher/=n_vec

    return fisher

def correlated_cl(par,mp) :
    nb=len(par.bins)-1
    lmax=3*par.nside-1
    pixsize=4*np.pi/hp.nside2npix(par.nside)
    cb=np.zeros(nb)
    larr=np.arange(lmax+1)

    z=invert_covar(par,mp)
    clz_scaled=(2*larr+1.)*hp.anafast(z)

    for ib in np.arange(nb) :
        cb[ib]=np.sum(clz_scaled[par.bins[ib]:par.bins[ib+1]])/pixsize**2

    return cb

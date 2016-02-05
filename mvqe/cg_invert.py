import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg,LinearOperator

glb_niter=0
def invert_covar(mp,ip_seen,signal_cl,sig_noise,ns,lmax_precon,tol,plot_stuff) :
    global glb_niter
    npx_seen=len(ip_seen)
    npx=hp.nside2npix(ns)
    pixsize=4*np.pi/npx
    glb_niter=0
    
    cl_signal_precon=pixsize**2*(1./signal_cl[:lmax_precon]-1./signal_cl[lmax_precon])
    sig_precon=np.sqrt(pixsize/signal_cl[lmax_precon])

    def linop_covar(m) :
        m_full=np.zeros(npx); m_full[ip_seen]=m
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0),signal_cl/pixsize),
                         ns,verbose=False)+sig_noise**2*m_full
        return m_out[ip_seen]

    def linop_precon(m) :
        m_full=np.zeros(npx); m_full[ip_seen]=m
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0,lmax=lmax_precon-1),
                                   cl_signal_precon),ns,verbose=False)/pixsize+sig_precon**2*m_full
        return m_out[ip_seen]

    def callback(m) :
        global glb_niter
        glb_niter+=1
        print glb_niter

    A=LinearOperator((npx_seen,npx_seen),matvec=linop_covar)
    M=LinearOperator((npx_seen,npx_seen),matvec=linop_precon)
    b=mp[ip_seen]

    m_solve,info=cg(A,b,x0=linop_precon(b),M=M,tol=tol,callback=callback)
    mp_out=np.zeros_like(mp); mp_out[ip_seen]=m_solve
    mp_fill=np.zeros_like(mp); mp_fill[ip_seen]=linop_covar(m_solve)

    plt.show()
    if plot_stuff :
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
        v_out[np.where(vr<0)]=-1
        return v_out

def fisher_single(v,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff,bins) :
    # v -> seen pixels
    nb=len(bins)-1
    npix=hp.nside2npix(nside)
    npix_seen=len(ip_seen)
    lmax=3*nside-1
    larr=np.arange(lmax+1)
    fisher=np.zeros([nb,nb])
    pixsize=4*np.pi/hp.nside2npix(nside)
    
    v_map=np.zeros(npix); v_map[ip_seen]=v
    vcm1=invert_covar(v_map,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff) #C^-1 * v
    v_lm=hp.map2alm(v_map,iter=0)
    vcm1_lm=hp.map2alm(vcm1,iter=0)
    for iba in np.arange(nb) :
        transfer=np.zeros(lmax+1); transfer[bins[iba]:bins[iba+1]]=1.
        v_map2=hp.alm2map(hp.almxfl(v_lm,transfer),nside,verbose=False)/pixsize #Q_a * v
        v_map2cm1=invert_covar(v_map2,ip_seen,signal_cl,sig_noise,
                               nside,lmax_precon,tol,plot_stuff) #C^-1 * Q_a * v
        va_lm=hp.map2alm(v_map2cm1,iter=0)
        cl_vcm1_va=(2*larr+1)*hp.alm2cl(vcm1_lm,alms2=va_lm)
        for ibb in np.arange(nb-iba)+iba :
            fisher[iba,ibb]=np.sum(cl_vcm1_va[bins[ibb]:bins[ibb+1]])/pixsize**2
            if iba!=ibb :
                fisher[ibb,iba]=fisher[iba,ibb]

    return fisher

def fisher_avg(n_vec,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff,bins) :
    nb=len(bins)-1
    npix_seen=len(ip_seen)
    fisher=np.zeros([nb,nb])
    
    for iv in np.arange(n_vec) :
        print "Realization #%d"%iv
        v=random_unit_variance(npix_seen,True)
        fisher+=fisher_single(v,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff,bins)
    fisher/=n_vec

    return fisher

def correlated_cl(mp,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff,bins) :
    nb=len(bins)-1
    lmax=3*nside-1
    pixsize=4*np.pi/hp.nside2npix(nside)
    cb=np.zeros(nb)
    larr=np.arange(lmax+1)

    z=invert_covar(mp,ip_seen,signal_cl,sig_noise,nside,lmax_precon,tol,plot_stuff)
    clz_scaled=(2*larr+1.)*hp.anafast(z)

    for ib in np.arange(nb) :
        cb[ib]=np.sum(clz_scaled[bins[ib]:bins[ib+1]])/pixsize**2

    return cb

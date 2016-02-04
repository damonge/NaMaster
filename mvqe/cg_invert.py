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
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0),signal_cl),
                         ns,verbose=False)/pixsize+sig_noise**2*m_full
        return m_out[ip_seen]

    def linop_precon(m) :
        m_full=np.zeros(npx); m_full[ip_seen]=m
        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0,lmax=lmax_precon-1),
                                   cl_signal_precon),ns,verbose=False)/pixsize+sig_precon**2*m_full
#        m_out=hp.alm2map(hp.almxfl(hp.map2alm(m_full,iter=0),
#                                   1./signal_cl),ns,verbose=False)*pixsize
        return m_out[ip_seen]

    def callback(m) :
        global glb_niter
        glb_niter+=1
        print glb_niter
#        m_fill=np.zeros_like(mp); m_fill[ip_seen]=linop_covar(m);
#        hp.mollview((mp-m_fill)/np.std(mp[ip_seen]))

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

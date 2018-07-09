from __future__ import print_function
from optparse import OptionParser
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pymaster as nmt
import os
import sys

DTOR=np.pi/180

def opt_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))
parser = OptionParser()
parser.add_option('--nside', dest='nside_out', default=512, type=int,
                  help='Resolution parameter')
parser.add_option('--isim-ini', dest='isim_ini', default=1, type=int,
                  help='Index of first simulation')
parser.add_option('--isim-end', dest='isim_end', default=100, type=int,
                  help='Index of last simulation')
parser.add_option('--wo-contaminants', dest='wo_cont', default=False, action='store_true',
                  help='Set if you don\'t want to use contaminants')
parser.add_option('--plot', dest='plot_stuff', default=False, action='store_true',
                  help='Set if you want to produce plots')
parser.add_option('--wo-pureb', dest='wo_pureb', default=False, action='store_true',
                  help='Set if you don\'t want to purify B-modes')
(o, args) = parser.parse_args()

nsims=o.isim_end-o.isim_ini+1
w_cont=not o.wo_cont
w_pureb=not o.wo_pureb

#Create output directory
predir="tests_sph"
os.system("mkdir -p "+predir)
prefix=predir+"/run_pure0%d_ns%d_cont%d"%(w_pureb,o.nside_out,w_cont)

#Read theory power spectra
def read_cl_camb(fname) :
    data=np.loadtxt(fname,unpack=True)
    ll=np.arange(3*o.nside_out,dtype=float)
    fac=2*np.pi/(ll[2:]*(ll[2:]+1.))
    cl_tt=np.zeros_like(ll); cl_tt[2:]=data[1,:3*o.nside_out-2]*fac
    cl_ee=np.zeros_like(ll); cl_ee[2:]=data[2,:3*o.nside_out-2]*fac
    cl_bb=np.zeros_like(ll); cl_bb[2:]=data[3,:3*o.nside_out-2]*fac
    cl_te=np.zeros_like(ll); cl_te[2:]=data[4,:3*o.nside_out-2]*fac

    return ll,cl_tt,cl_ee,cl_bb,cl_te
l,cltt,clee,clbb,clte=read_cl_camb("data/cls_cmb.txt")
#Noise power spectrum
nlev=(1.*np.pi/(180*60))**2 #1 uK-arcmin noise level
nltt=nlev*(np.ones_like(l)+(10./(l+0.1))**2.4) #1/ell noise with a knee scale of ell=10 (optimistic)
nlee=2*nltt; nlbb=2*nltt; nlte=0*nltt
#Beam
fwhm_amin=25. #Corresponding to 0.5m aperture at 90GHz
beam=np.exp(-0.5*l*(l+1)*(fwhm_amin*np.pi/(180*60*2.355))**2)
if o.plot_stuff :
    plt.figure()
    plt.plot(l,clee,'r-',label='EE signal')
    plt.plot(l,clbb,'b-',label='BB signal')
    plt.plot(l,nlee/beam**2,'k--',label='Noise')
    plt.loglog(); plt.legend()
    plt.xlabel('$\\ell$',fontsize=16);
    plt.ylabel('$C_\\ell$',fontsize=16);

#Read mask
mask=hp.read_map("data/mask_cmb_ns%d.fits"%o.nside_out,verbose=False)
if o.plot_stuff :
    hp.mollview(mask)
fsky=np.mean(mask/np.amax(mask));

#Read contaminant maps
if w_cont :
    fgp=np.zeros([1,2,hp.nside2npix(o.nside_out)]);
    fgp[0,0,:],fgp[0,1,:]=hp.read_map("data/cont_cmb_ns%d.fits"%o.nside_out,
                                      field=[0,1],verbose=False); #Foregrounds
    fgp[0,0,:]=hp.smoothing(fgp[0,0,:],beam=beam,verbose=False)
    fgp[0,1,:]=hp.smoothing(fgp[0,1,:],beam=beam,verbose=False)
    if o.plot_stuff :
        hp.mollview(np.sum(fgp,axis=0)[0,:]*mask)
        hp.mollview(np.sum(fgp,axis=0)[1,:]*mask)

#Binning scheme
ls=np.arange(3*o.nside_out,dtype=int)
bpws=np.zeros(3*o.nside_out,dtype=int)-1
weights=np.ones(3*o.nside_out)
bpw_edges=[2,9,17]
while bpw_edges[-1]<3*o.nside_out :
    bpw_edges.append(min(bpw_edges[-1]+12,3*o.nside_out))
bpw_edges=np.array(bpw_edges)
for ib,b0 in enumerate(bpw_edges[:-1]) :
    bpws[b0:bpw_edges[ib+1]]=ib
    weights[b0:bpw_edges[ib+1]]=1./(bpw_edges[ib+1]-b0+0.)
b=nmt.NmtBin(o.nside_out,ells=ls,bpws=bpws,weights=weights)

#Generate some initial fields
print(" - Res: %.3lf arcmin. "%(np.sqrt(4*np.pi*(180*60/np.pi)**2/hp.nside2npix(o.nside_out))))
def get_fields() :
    #Signal
    st,sq,su=hp.synfast([cltt*beam**2+nltt,clee*beam**2+nlee,
                         clbb*beam**2+nlbb,clte*beam**2+nlte],
                        o.nside_out,new=True,verbose=False,pol=True)
    if w_cont :
        sq+=np.sum(fgp,axis=0)[0,:]; su+=np.sum(fgp,axis=0)[1,:]
        ff2=nmt.NmtField(mask,[sq,su],templates=fgp,beam=beam,
                         purify_e=False,purify_b=w_pureb,n_iter_mask_purify=10)
    else :
        ff2=nmt.NmtField(mask,[sq,su],beam=beam,
                         purify_e=False,purify_b=w_pureb,n_iter_mask_purify=10)
    return ff2
np.random.seed(1000)
f2=get_fields()
    
if o.plot_stuff :
    hp.mollview(f2.get_maps()[0]*mask,title='$Q$')
    hp.mollview(f2.get_maps()[1]*mask,title='$U$')

#Use initial fields to generate coupling matrix
w22=nmt.NmtWorkspace();
if not os.path.isfile(prefix+"_w22.dat") :
    print("Computing 22")
    w22.compute_coupling_matrix(f2,f2,b)
    w22.write_to(prefix+"_w22.dat");
else :
    w22.read_from(prefix+"_w22.dat")

#Generate theory prediction
cl22_th=w22.decouple_cell(w22.couple_cell([clee,0*clee,0*clbb,clbb]))
np.savetxt(prefix+"_cl_th.txt",
           np.transpose([b.get_effective_ells(),cl22_th[0],cl22_th[1],cl22_th[2],cl22_th[3]]))

#Compute noise and deprojection bias
if not os.path.isfile(prefix+"_clb22.npy") :
    print("Computing deprojection and noise bias 00")
    #Compute noise bias
    clb22=w22.couple_cell([nlee,0*nlee,0*nlbb,nlbb])
    #Compute deprojection bias
    if w_cont :
        #Signal contribution
        clb22+=nmt.deprojection_bias(f0,f0,[clee*beam**2+nlee,0*clee,0*clbb,clbb*beam**2+nlbb])
    np.save(prefix+"_clb22",clb22)
else :
    clb22=np.load(prefix+"_clb22.npy")
    
#Compute mean and variance over nsims simulations
cl22_all=[]
for i in np.arange(nsims) :
    #if i%100==0 :
    print("%d-th sim"%(i+o.isim_ini))

    if not os.path.isfile(prefix+"_cl_%04d.txt"%(o.isim_ini+i)) :
        np.random.seed(1000+o.isim_ini+i)
        f2=get_fields()
        cl22=w22.decouple_cell(nmt.compute_coupled_cell(f2,f2),cl_bias=clb22)
        np.savetxt(prefix+"_cl_%04d.txt"%(o.isim_ini+i),
                   np.transpose([b.get_effective_ells(),cl22[0],cl22[1],cl22[2],cl22[3]]))
    cld=np.loadtxt(prefix+"_cl_%04d.txt"%(o.isim_ini+i),unpack=True)
    cl22_all.append([cld[1],cld[2],cld[3],cld[4]])
cl22_all=np.array(cl22_all)

#Plot results
if o.plot_stuff :
    l_eff=b.get_effective_ells()
    cols=plt.cm.rainbow(np.linspace(0,1,3))
    plt.figure()
    plt.errorbar(l_eff,np.mean(cl22_all,axis=0)[0]/cl22_th[0]-1,
                 yerr=np.std(cl22_all,axis=0)[0]/cl22_th[0]/np.sqrt(nsims+0.),
                 label='$EE$',fmt='bo')
    plt.errorbar(l_eff,np.mean(cl22_all,axis=0)[3]/cl22_th[3]-1,
                 yerr=np.std(cl22_all,axis=0)[3]/cl22_th[3]/np.sqrt(nsims+0.),
                 label='$BB$',fmt='ro')
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$\\Delta C_\\ell/C_\\ell$',fontsize=16)
    plt.xlim([2,2*o.nside_out])
    plt.ylim([-0.2,0.2])
    plt.legend(loc='lower right',frameon=False,fontsize=16)
    plt.savefig(prefix+'_celldiff.png',bbox_inches='tight')
    plt.savefig(prefix+'_celldiff.pdf',bbox_inches='tight')

    import scipy.stats as st
    bins_use=np.where(l_eff<2*o.nside_out)[0]; ndof=len(bins_use)
    res=(cl22_all[:,:,:]-cl22_th[None,:,:])/np.std(cl22_all,axis=0)[None,:,:]
    chi2_22=np.sum(res[:,:,bins_use]**2,axis=2)

    x=np.linspace(ndof-5*np.sqrt(2.*ndof),ndof+5*np.sqrt(2*ndof),256)
    pdf=st.chi2.pdf(x,ndof)

    plt.figure(figsize=(10,4))
    ax=[plt.subplot(1,3,i+1) for i in range(3)]
    plt.subplots_adjust(wspace=0, hspace=0)

    h,b,p=ax[0].hist(chi2_22[:,0],bins=40,density=True)
    ax[0].text(0.75,0.9,'$EE$',transform=ax[0].transAxes)
    ax[0].set_xlabel('$\\chi^2$')
    ax[0].set_ylabel('$P(\\chi^2)$')

    h,b,p=ax[1].hist(chi2_22[:,1],bins=40,density=True)
    ax[1].text(0.75,0.9,'$EB$',transform=ax[1].transAxes)

    h,b,p=ax[2].hist(chi2_22[:,3],bins=40,density=True)
    ax[2].text(0.75,0.9,'$BB$',transform=ax[2].transAxes)

    for a in ax :
        a.set_xlabel('$\\chi^2$')
    ax[1].set_yticklabels([])
    ax[2].set_yticklabels([])
    for a in ax :
        a.set_xlim([ndof-5*np.sqrt(2.*ndof),ndof+5*np.sqrt(2.*ndof)])
        a.set_ylim([0,1.4*np.amax(pdf)])
        a.plot([ndof,ndof],[0,1.4*np.amax(pdf)],'k--',label='$N_{\\rm dof}$')
        a.plot(x,pdf,'k-',label='$P(\\chi^2,N_{\\rm dof})$')
    ax[0].legend(loc='upper left',frameon=False)
    plt.savefig(prefix+'_distributions.png',bbox_inches='tight')
    plt.savefig(prefix+'_distributions.pdf',bbox_inches='tight')

    ic=0
    plt.figure()
    plt.plot(l_eff,np.mean(cl22_all,axis=0)[0],'.',
             label='$EE$',c=cols[ic]);
    plt.plot(l_eff,cl22_th[0],'--',c=cols[ic]); ic+=1
    plt.plot(l_eff,np.mean(cl22_all,axis=0)[1],'.',
             label='$EB$',c=cols[ic]); ic+=1
    plt.plot(l_eff,np.mean(cl22_all,axis=0)[3],'.',
             label='$BB$',c=cols[ic]);
    plt.plot(l_eff,cl22_th[3],'--',c=cols[ic]); ic+=1
    plt.yscale('log')
    plt.xlim([2,2*o.nside_out])
    plt.xlabel('$\\ell$',fontsize=16)
    plt.ylabel('$C_\\ell$',fontsize=16)
    plt.legend(loc='lower left',frameon=False,fontsize=14,ncol=2)
    plt.savefig(prefix+'_cellfull.png',bbox_inches='tight')
    plt.savefig(prefix+'_cellfull.pdf',bbox_inches='tight')
    plt.show()

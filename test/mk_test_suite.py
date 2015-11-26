import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

la,dl_ttr,dl_eer,dl_bbr,dl_ter=np.loadtxt("cl.txt",unpack=True)
dl_tt=np.zeros(len(la)+2); dl_tt[2:]=dl_ttr
dl_ee=np.zeros(len(la)+2); dl_ee[2:]=dl_eer
dl_bb=np.zeros(len(la)+2); dl_bb[2:]=dl_bbr
dl_te=np.zeros(len(la)+2); dl_te[2:]=dl_ter
l=np.arange(len(dl_tt))

cl_tt=dl_tt*2*np.pi/(l*(l+1)); cl_tt[0]=0;
cl_ee=dl_ee*2*np.pi/(l*(l+1)); cl_ee[0]=0;
cl_bb=dl_bb*2*np.pi/(l*(l+1)); cl_bb[0]=0;
cl_te=dl_te*2*np.pi/(l*(l+1)); cl_te[0]=0;
cl_tb=np.zeros_like(cl_te)
cl_eb=np.zeros_like(cl_te)
dl_tb=np.zeros_like(cl_te)
dl_eb=np.zeros_like(cl_te)

nside=256

alms=hp.synalm([cl_tt,cl_ee,cl_bb,cl_te,cl_eb,cl_tb],lmax=3*nside-1,new=True)
map_t,map_e,map_b=hp.alm2map(alms,nside,pol=False)
map_tb,map_q,map_u=hp.alm2map(alms,nside,pol=True)

#hp.mollview(map_t);
#hp.mollview(map_e);
#hp.mollview(map_b);
#hp.mollview(map_tb);
#hp.mollview(map_q);
#hp.mollview(map_u);
#plt.show()

cl_tt_a,cl_ee_a,cl_bb_a,cl_te_a,cl_eb_a,cl_tb_a=hp.anafast([map_t,map_e,map_b],pol=False)
cl_tt_b,cl_ee_b,cl_bb_b,cl_te_b,cl_eb_b,cl_tb_b=hp.anafast([map_tb,map_q,map_u],pol=True)
l_d=np.arange(len(cl_tt_a));
plt.plot(l[:len(l_d)],dl_tt[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_tt_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_tt_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.gca().set_yscale('log'); plt.show()
plt.plot(l[:len(l_d)],dl_ee[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_ee_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_ee_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.gca().set_yscale('log'); plt.show()
plt.plot(l[:len(l_d)],dl_bb[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_bb_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_bb_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.gca().set_yscale('log'); plt.show()

plt.plot(l[:len(l_d)],dl_te[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_te_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_te_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.show()
plt.plot(l[:len(l_d)],dl_eb[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_eb_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_eb_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.show()
plt.plot(l[:len(l_d)],dl_tb[:len(l_d)],'k-'); plt.plot(l_d,l_d*(l_d+1)*cl_tb_a/(2*np.pi),'r-'); plt.plot(l_d,l_d*(l_d+1)*cl_tb_b/(2*np.pi),'b-'); plt.gca().set_xscale('log'); plt.show()

mask=np.ones(hp.nside2npix(nside))
hp.write_map("mask_full.fits",mask);
th,phi=hp.pix2ang(nside,np.arange(hp.nside2npix(nside)))
dec=90-180*th/np.pi
mask[np.where(np.fabs(dec)<10.)]=0
hp.write_map("mask_cut.fits",mask);

hp.write_map("map_t.fits",map_t,column_names=["t"])
hp.write_map("map_e.fits",map_e,column_names=["e"])
hp.write_map("map_b.fits",map_b,column_names=["b"])
hp.write_map("map_qu.fits",[map_q,map_u],column_names=["q","u"]);
noise=np.zeros_like(l_d)
np.savetxt("cl_noise_1.dat",np.transpose([l_d,noise]))
np.savetxt("cl_noise_2.dat",np.transpose([l_d,noise,noise]))
np.savetxt("cl_noise_4.dat",np.transpose([l_d,noise,noise,noise,noise]))

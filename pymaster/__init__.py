"""
:mod:`pymaster` contains three basic classes:

- :class:`pymaster.NmtField`
- :class:`pymaster.NmtBin`
- :class:`pymaster.NmtWorkspace`

And a number of functions

- :func:`pymaster.deprojection_bias`
- :func:`pymaster.compute_coupled_cell`
- :func:`pymaster.compute_full_master`
- :func:`pymaster.mask_apodization`

Many of these function accept or return sets of power spectra (arrays with one element per angular multipole) or bandpowers (binned versions of power spectra). In all cases, these are returned and provided as 2D arrays with shape [n_cls][nl], where n_cls is the number of power spectra and nl is either the number of multipoles or bandpowers. In all cases, n_cls should correspond with the spins of the two fields being correlated, and the ordering is as follows:

- Two spin-0 fields: n_cls=1, [C_T1T2]
- One spin-0 field and one spin-2 field: n_cls=2, [C_TE,C_TB]
- Two spin-2 fields: n_cls=4, [C_E1E2,C_E1B2,C_E2B1,C_B1B2]

All sky maps accepted and returned by these functions are in the form of HEALPix maps exclusively with RING ordering.

"""
import nmtlib as lib
import numpy as np


class NmtField(object) :
    """
    An NmtField object contains all the information describing the fields to correlate, including their observed maps, masks and contaminant templates.

    :param mask: array containing a HEALPix map corresponding to the field's mask.
    :param maps: 2D array containing the observed maps for this field. The first dimension corresponds to the number of maps, which should be 1 for a spin-0 field and 2 for a spin-2 field.
    :param templates: 3D array containing a set of contaminant templates for this field. This array should have shape [ntemp][nmap][npix], where ntemp is the number of templates, nmap should be 1 for spin-0 fields and 2 for spin-2 fields, and npix is the number of pixels per map. The best-fit contribution from each contaminant is automatically removed from the maps unless templates=None
    :param beam: spherical harmonic transform of the instrumental beam (assumed to be rotationally symmetric - i.e. no m dependence). If None, no beam will be corrected for. Otherwise, this array should have 3*nside elements, corresponding to multipoles from 0 to 3*nside-1.
    :param purify_e: use pure E-modes?
    :param purify_b: use pure B-modes?

    """
    def __init__(self,mask,maps,templates=None,beam=None,purify_e=False,purify_b=False,lx=-1.,ly=-1.,lmax=-1) :
        pure_e=0
        if(purify_e) :
            pure_e=1
        pure_b=0
        if(purify_b) :
            pure_b=1

        self.isflat=False
        if len(mask)==2 :
            self.isflat=True

        if self.isflat :
            if (lx<0) or (ly<0) :
                raise KeyError("Must supply dimensions for flat-sky field")
            #Flatten arrays and check dimensions
            shape_2D=np.shape(mask)
            nmaps=len(maps[0])
            self.ny=shape_2D[0]
            self.nx=shape_2D[1]

            #Flatten mask
            msk=mask.flatten()

            #Flatten maps
            mps=[]
            for m in maps :
                if np.shape(m)!=shape_2d :
                    KeyError("Mask and maps don't have the same shape")
                mps.append(m.flatten())
            mps=np.array(mps)

            #Flatten templates
            if templates!=None :
                tmps=[]
                for t in templates :
                    tmp=[]
                    if len(t)!=nmaps :
                        KeyError("Maps and templates should have the same number of maps")
                    for m in t :
                        tmp.append(m.flatten())
                    tmps.append(tmp)
                tmps=np.array(tmps)

            #Form beam
            if(beam==None) :
                if lmax<=0 :
                    KeyError("Lmax must be positive")
                beam_use=np.ones(lmax+1)
            else :
                beam_use=beam
                lmax=len(beam_use)-1

            #Generate field
            if(templates==None) :
                self.fl=lib.field_alloc_new_notemp_flat(nx,ny,lx,ly,msk,mps,beam_use,pure_e,pure_b)
            else :
                self.fl=lib.field_alloc_new_flat(nx,ny,lx,ly,msk,mps,tmps,beam_use,pure_e,pure_b)
        else :
            nside=2
            while(12*nside*nside!=len(mask)) :
                nside*=2
                if(nside>65536) :
                    raise KeyError("Something is wrong with your input arrays")

            if((len(maps)!=1) and (len(maps)!=2)) :
                raise KeyError("Must supply 1 or 2 maps per field")
            if(templates!=None) :
                if((len(templates[0])!=1) and (len(templates[0])!=2)) :
                    raise KeyError("Must supply 1 or 2 maps per field")
            if(beam==None) :
                beam_use=np.ones(3*nside)
            else :
                if(len(beam)!=3*nside) :
                    raise KeyError("Input beam must have 3*nside elements")
                beam_use=beam

            if(templates==None) :
                self.fl=lib.field_alloc_new_notemp(mask,maps,beam_use,pure_e,pure_b)
            else :
                self.fl=lib.field_alloc_new(mask,maps,templates,beam_use,pure_e,pure_b)

    def __del__(self) :
        lib.field_free(self.fl)

    def get_maps(self) :
        """
        Returns a 2D array ([nmap][npix]) corresponding to the observed maps for this field. If the field was initialized with contaminant templates, the maps returned by this function have their best-fit contribution from these contaminants removed.
        
        :return: 2D array of HEALPix maps
        """
        maps=np.zeros([self.fl.nmaps,self.fl.npix])
        for imap in np.arange(self.fl.nmaps) :
            maps[imap,:]=lib.get_map(self.fl,imap,int(self.fl.npix))
        if self.isflat :
            mps=maps.reshape([self.fl.nmaps,self.ny,self.nx])
        else :
            mps=maps
        return mps

    def get_templates(self) :
        """
        Returns a 3D array ([ntemp][nmap][npix]) corresponding to the contaminant templates passed when initializing this field.

        :return: 3D array of HEALPix maps
        """
        temp=np.zeros([self.fl.ntemp,self.fl.nmaps,self.fl.npix])
        for itemp in np.arange(self.fl.ntemp) :
            for imap in np.arange(self.fl.nmaps) :
                temp[itemp,imap,:]=lib.get_temp(self.fl,itemp,imap,int(self.fl.npix))
        if self.isflat :
            tmps=temp.reshape([self.fl.ntemp,self.fl.nmaps,self.ny,self.nx])
        else :
            tmps=temp
        return tmps


class NmtBin(object) :
    """
    An NmtBin object defines the set of bandpowers used in the computation of the pseudo-Cl estimator. The definition of bandpowers is described in Section 3.6 of the scientific documentation.

    :param int nside: HEALPix nside resolution parameter of the maps you intend to correlate. The maximum multipole considered for bandpowers will be 3*nside-1.
    :param array-like ells: array of integers corresponding to different multipoles
    :param array-like bpws: array of integers that assign the multipoles in ells to different bandpowers
    :param array-like weights: array of floats corresponding to the weights associated to each multipole in ells. The sum of weights within each bandpower is normalized to 1.
    :param int nlb: integer value corresponding to a constant bandpower width. I.e. the bandpowers will be defined as consecutive sets of nlb multipoles from l=2 to l=3*nside-1 with equal weights. If this argument is provided, the values of ells, bpws and weights are ignored.
    """
    def __init__(self,nside,bpws=None,ells=None,weights=None,nlb=None) :
        if((bpws==None) and (ells==None) and (weights==None) and (nlb==None)) :
            raise KeyError("Must supply bandpower arrays or constant bandpower width")

        if(nlb==None) :
            if((bpws==None) and (ells==None) and (weights==None)) :
                raise KeyError("Must provide bpws, ells and weights")
            self.bin=lib.bins_create_py(bpws,ells,weights,3*nside-1)
        else :
            self.bin=lib.bins_constant(nlb,3*nside-1)
        self.lmax=3*nside-1

    def __del__(self) :
        lib.bins_free(self.bin)

    def get_n_bands(self) :
        """
        Returns the number of bandpowers stored in this object

        :return: number of bandpowers
        """
        return self.bin.n_bands

    def get_nell_list(self) :
        """
        Returns an array with the number of multipoles in each bandpower stored in this object

        :return: number of multipoles per bandpower
        """
        return lib.get_nell_list(self.bin,self.bin.n_bands)

    def get_ell_list(self,ibin) :
        """
        Returns an array with the multipoles in the ibin-th bandpower

        :param int ibin: bandpower index
        :return: multipoles associated with bandpower ibin
        """
        return lib.get_ell_list(self.bin,ibin,lib.get_nell(self.bin,ibin))

    def get_weight_list(self,ibin) :
        """
        Returns an array with the weights associated to each multipole in the ibin-th bandpower

        :param int ibin: bandpower index
        :return: weights associated to multipoles in bandpower ibin
        """
        return lib.get_weight_list(self.bin,ibin,lib.get_nell(self.bin,ibin))

    def get_effective_ells(self) :
        """
        Returns an array with the effective multipole associated to each bandpower. These are computed as a weighted average of the multipoles within each bandpower.

        :return: effective multipoles for each bandpower
        """
        return lib.get_ell_eff(self.bin,self.bin.n_bands)

    def bin_cell(self,cls_in) :
        """
        Bins a power spectrum into bandpowers. This is carried out as a weighted average over the multipoles in each bandpower.

        :param array-like cls_in: 2D array of power spectra
        :return: array of bandpowers
        """
        if(len(cls_in[0])!=self.lmax+1) :
            raise KeyError("Input Cl has wrong size")
        cl1d=lib.bin_cl(self.bin,cls_in,len(cls_in)*self.bin.n_bands)
        clout=np.reshape(cl1d,[len(cls_in),self.bin.n_bands])
        return clout

    def unbin_cell(self,cls_in) :
        """
        Un-bins a set of bandpowers into a power spectrum. This is simply done by assigning a constant value for every multipole in each bandpower (corresponding to the value of that bandpower).

        :param array-like cls_in: array of bandpowers
        :return: array of power spectra
        """
        if(len(cls_in[0])!=self.bin.n_bands) :
            raise KeyError("Input Cl has wrong size")
        cl1d=lib.unbin_cl(self.bin,cls_in,len(cls_in)*(self.lmax+1))
        clout=np.reshape(cl1d,[len(cls_in),self.lmax+1])
        return clout


class NmtWorkspace(object) :
    """
    NmtWorkspace objects are used to compute and store the coupling matrix associated with an incomplete sky coverage, and used in the MASTER algorithm. When initialized, this object is practically empty. The information describing the coupling matrix must be computed or read from a file afterwards.
    """
    def __init__(self) :
        self.wsp=None

    def __del__(self) :
        if(self.wsp!=None) :
            lib.workspace_free(self.wsp)

    def read_from(self,fname) :
        """
        Reads the contents of an NmtWorkspace object from a file (encoded using an internal binary format).

        :param str fname: input file name
        """
        if self.wsp!=None :
            lib.workspace_free(self.wsp)
        self.wsp=lib.workspace_read(fname);
        
    def compute_coupling_matrix(self,fl1,fl2,bins) :
        """
        Computes coupling matrix associated with the cross-power spectrum of two NmtFields and an NmtBin binning scheme.

        :param NmtField fl1,fl2: fields to correlate
        :param NmtBin bin: binning scheme
        """
        if self.wsp!=None :
            lib.workspace_free(self.wsp)
        if(fl1.isflat and fl2.isflat) :
            raise KeyError("MASTER for flat-sky fields still not supported")
        self.wsp=lib.compute_coupling_matrix(fl1.fl,fl2.fl,bins.bin)

    def write_to(self,fname) :
        """
        Writes the contents of an NmtWorkspace object to a file (encoded using an internal binary format).

        :param str fname: output file name
        """
        if self.wsp==None :
            raise KeyError("Must initialize workspace before writing")
        lib.workspace_write(self.wsp,fname)

    def couple_cell(self,cl_in) :
        """
        Convolves a set of input power spectra with a coupling matrix (see Eq. 6 of the C API documentation).

        :param cl_in: set of input power spectra. The number of power spectra must correspond to the spins of the two fields that this NmtWorkspace object was initialized with (i.e. 1 for two spin-0 fields, 2 for one spin-0 and one spin-2 field and 4 for two spin-2 fields).
        :return: coupled power spectrum
        """
        if((len(cl_in)!=self.wsp.ncls) or (len(cl_in[0])!=self.wsp.lmax+1)) :
            raise KeyError("Input power spectrum has wrong shape")
        cl1d=lib.couple_cell_py(self.wsp,cl_in,self.wsp.ncls*(self.wsp.lmax+1))
        clout=np.reshape(cl1d,[self.wsp.ncls,self.wsp.lmax+1])
        return clout

    def decouple_cell(self,cl_in,cl_bias=None,cl_noise=None) :
        """
        Decouples a set of pseudo-Cl power spectra into a set of bandpowers by inverting the binned coupling matrix (se Eq. 4 of the C API documentation).

        :param cl_in: set of input power spectra. The number of power spectra must correspond to the spins of the two fields that this NmtWorkspace object was initialized with (i.e. 1 for two spin-0 fields, 2 for one spin-0 and one spin-2 field and 4 for two spin-2 fields).
        :param cl_bias: bias to the power spectrum associated to contaminant residuals (optional). This can be computed through :func:`pymaster.deprojection_bias`.
        :param cl_noise: noise bias (i.e. angular power spectrum of masked noise realizations).
        :return: set of decoupled bandpowers
        """
        if((len(cl_in)!=self.wsp.ncls) or (len(cl_in[0])!=self.wsp.lmax+1)) :
            raise KeyError("Input power spectrum has wrong shape")
        if cl_bias!=None :
            if((len(cl_bias)!=self.wsp.ncls) or (len(cl_bias[0])!=self.wsp.lmax+1)) :
                raise KeyError("Input bias power spectrum has wrong shape")
            clb=cl_bias.copy()
        else :
            clb=np.zeros_like(cl_in)
        if cl_noise!=None :
            if((len(cl_noise)!=self.wsp.ncls) or (len(cl_noise[0])!=self.wsp.lmax+1)) :
                raise KeyError("Input noise power spectrum has wrong shape")
            cln=cl_noise.copy()
        else :
            cln=np.zeros_like(cl_in)

        cl1d=lib.decouple_cell_py(self.wsp,cl_in,cln,clb,self.wsp.ncls*self.wsp.bin.n_bands)
        clout=np.reshape(cl1d,[self.wsp.ncls,self.wsp.bin.n_bands])

        return clout

def deprojection_bias(f1,f2,cls_guess) :
    """
    Computes the bias associated to contaminant removal to the cross-pseudo-Cl of two fields.

    :param NmtField f1,f2: fields to correlate
    :param cls_guess: set of power spectra corresponding to a best-guess of the true power spectra of f1 and f2.
    :return: deprojection bias power spectra.
    """
    if(fl1.isflat and fl2.isflat) :
        raise KeyError("MASTER for flat-sky fields still not supported")
    if(len(cls_guess)!=f1.fl.nmaps*f2.fl.nmaps) :
        raise KeyError("Proposal Cell doesn't match number of maps")
    if(len(cls_guess[0])!=f1.fl.lmax+1) :
        raise KeyError("Proposal Cell doesn't match map resolution")
    cl1d=lib.comp_deproj_bias(f1.fl,f2.fl,cls_guess,len(cls_guess)*len(cls_guess[0]))
    cl2d=np.reshape(cl1d,[len(cls_guess),len(cls_guess[0])])

    return cl2d

def compute_coupled_cell(f1,f2,n_iter=3) :
    """
    Computes the full-sky angular power spectra of two masked fields (f1 and f2) without aiming to deconvolve the mode-coupling matrix. Effectively, this is equivalent to calling the usual HEALPix anafast routine on the masked and contaminant-cleaned maps.

    :param NmtField f1,f2: fields to correlate
    :param int n_iter: number of iterations for SHTs (optional)
    :return: array of coupled power spectra
    """
    if(fl1.isflat and fl2.isflat) :
        raise KeyError("MASTER for flat-sky fields still not supported")
    if(f1.fl.nside!=f2.fl.nside) :
        raise KeyError("Fields must have same resolution")
    
    cl1d=lib.comp_pspec_coupled(f1.fl,f2.fl,f1.fl.nmaps*f2.fl.nmaps*(f1.fl.lmax+1),n_iter)
    clout=np.reshape(cl1d,[f1.fl.nmaps*f2.fl.nmaps,f1.fl.lmax+1])

    return clout
    
def compute_full_master(f1,f2,b,cl_noise=None,cl_guess=None,workspace=None) :
    """
    Computes the full MASTER estimate of the power spectrum of two fields (f1 and f2). This is equivalent to successively calling:

    - :func:`pymaster.NmtWorkspace.compute_coupling_matrix`
    - :func:`pymaster.deprojection_bias`
    - :func:`pymaster.compute_coupled_cell`
    - :func:`pymaster.NmtWorkspace.decouple_cell`

    :param NmtField f1,f2: fields to correlate
    :param NmtBin b: binning scheme defining output bandpower
    :param cl_noise: noise bias (i.e. angular power spectrum of masked noise realizations) (optional).
    :param cls_guess: set of power spectra corresponding to a best-guess of the true power spectra of f1 and f2. Needed only to compute the contaminant cleaning bias (optional).
    :param NmtWorkspace workspace: object containing the mode-coupling matrix associated with an incomplete sky coverage. If provided, the function will skip the computation of the mode-coupling matrix and use the information encoded in this object.
    :return: set of decoupled bandpowers
    """
    if(fl1.isflat and fl2.isflat) :
        raise KeyError("MASTER for flat-sky fields still not supported")
    if(f1.fl.nside!=f2.fl.nside) :
        raise KeyError("Fields must have same resolution")
    if cl_noise!=None :
        if(len(cl_noise)!=f1.fl.nmaps*f2.fl.nmaps) :
            raise KeyError("Wrong length for noise power spectrum")
        cln=cl_noise.copy()
    else :
        cln=np.zeros([f1.fl.nmaps*f2.fl.nmaps,3*f1.fl.nside])
    if cl_guess!=None :
        if(len(cl_guess)!=f1.fl.nmaps*f2.fl.nmaps) :
            raise KeyError("Wrong length for guess power spectrum")
        clg=cl_guess.copy()
    else :
        clg=np.zeros([f1.fl.nmaps*f2.fl.nmaps,3*f1.fl.nside])

    if workspace==None :
        cl1d=lib.comp_pspec(f1.fl,f2.fl,b.bin,None,cln,clg,len(cln)*b.bin.n_bands)
    else :
        cl1d=lib.comp_pspec(f1.fl,f2.fl,b.bin,workspace.wsp,cln,clg,len(cln)*b.bin.n_bands)

    clout=np.reshape(cl1d,[len(cln),b.bin.n_bands])

    return clout

def mask_apodization(mask_in,aposize,apotype="C1") :
    """
    Apodizes a mask with an given apodization scale using different methods.

    :param mask_in: input mask, provided as an array of plots corresponding to a HEALPix map in RING order.
    :param aposize: apodization scale in degrees.
    :param apotype: apodization type. Three methods implemented: "C1", "C2" and "Smooth". See the description of the C-function nmt_apodize_mask in the C API documentation for a full description of these methods.
    :return: apodized mask as a HEALPix map
    """
    return lib.apomask(mask_in,len(mask_in),aposize,apotype)

def synfast_spherical(nside,cls,pol=False,beam=None) :
    seed=np.random.randint(50000000)
    if pol :
        use_pol=1
        nmaps=3
        if(len(np.shape(cls))!=2) :
            raise KeyError("You should supply more than one power spectrum if you want polarization")
        ncl=len(cls)
        if ((ncl!=4) and (ncl!=6)) :
            raise KeyError("You should provide 4 or 6 power spectra if you want polarization")
        lmax=len(cls[0])-1
        cls_use=np.zeros([6,lmax+1])
        cls_use[0,:]=cls[0] #TT
        cls_use[1,:]=cls[3] #TE
        cls_use[3,:]=cls[1] #EE
        cls_use[5,:]=cls[2] #BB
        if ncl==6 :
            cls_use[2,:]=cls[4] #TB
            cls_use[4,:]=cls[5] #EB
    else :
        use_pol=0
        nmaps=1
        if(len(np.shape(cls))!=1) :
            raise KeyError("You should supply only one power spectrum if you don't want polarization")
        lmax=len(cls)-1
        cls_use=np.array([cls])

    if beam==None :
        beam_use=np.ones(lmax+1)
    else :
        if len(beam)!=lmax+1 :
            raise KeyError("The beam should have as many multipoles as the power spectrum")
        beam_use=beam
    data=lib.synfast_new(nside,use_pol,seed,cls_use,beam_use,nmaps*12*nside*nside)

    maps=data.reshape([nmaps,12*nside*nside])

    return maps

def synfast_flat(nx,ny,lx,ly,cls,pol=False,beam=None) :
    seed=np.random.randint(50000000)
    if pol :
        use_pol=1
        nmaps=3
        if(len(np.shape(cls))!=2) :
            raise KeyError("You should supply more than one power spectrum if you want polarization")
        ncl=len(cls)
        if ((ncl!=4) and (ncl!=6)) :
            raise KeyError("You should provide 4 or 6 power spectra if you want polarization")
        lmax=len(cls[0])-1
        cls_use=np.zeros([6,lmax+1])
        cls_use[0,:]=cls[0] #TT
        cls_use[1,:]=cls[3] #TE
        cls_use[3,:]=cls[1] #EE
        cls_use[5,:]=cls[2] #BB
        if ncl==6 :
            cls_use[2,:]=cls[4] #TB
            cls_use[4,:]=cls[5] #EB
    else :
        use_pol=0
        nmaps=1
        if(len(np.shape(cls))!=1) :
            raise KeyError("You should supply only one power spectrum if you don't want polarization")
        lmax=len(cls)-1
        cls_use=np.array([cls])

    if beam==None :
        beam_use=np.ones(lmax+1)
    else :
        if len(beam)!=lmax+1 :
            raise KeyError("The beam should have as many multipoles as the power spectrum")
        beam_use=beam
    data=lib.synfast_new_flat(nx,ny,lx,ly,use_pol,seed,cls_use,beam_use,nmaps*ny*nx)

    maps=data.reshape([nmaps,ny,nx])

    return maps

import nmtlib as lib
import numpy as np

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
    if(f1.isflat and f2.isflat) :
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
    if(f1.isflat and f2.isflat) :
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
    if(f1.isflat and f2.isflat) :
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

import nmtlib as lib
import numpy as np

class NmtWorkspace(object) :
    def __init__(self) :
        self.wsp=None

    def __del__(self) :
        if(self.wsp!=None) :
            lib.workspace_free(self.wsp)

    def read_from(self,fname) :
        if self.wsp!=None :
            lib.workspace_free(self.wsp)
        self.wsp=lib.workspace_read(fname);
        
    def compute_coupling_matrix(self,fl1,fl2,bins) :
        if self.wsp!=None :
            lib.workspace_free(self.wsp)
        self.wsp=lib.compute_coupling_matrix(fl1.fl,fl2.fl,bins.bin)

    def write_to(self,fname) :
        if self.wsp==None :
            raise KeyError("Must initialize workspace before writing")
        lib.workspace_write(self.wsp,fname)

    def couple(self,cl_in) :
        if((len(cl_in)!=self.wsp.ncls) or (len(cl_in[0])!=self.wsp.lmax+1)) :
            raise KeyError("Input power spectrum has wrong shape")
        cl1d=lib.couple_cell(self.wsp,cl_in,self.wsp.ncls*(self.wsp.lmax+1))
        clout=np.reshape(cl1d,[self.wsp.ncls,self.wsp.lmax+1])
        return clout

    def decouple_cell(self,cl_in,cl_bias=None,cl_noise=None) :
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

        cl1d=lib.decouple_cell(self.wsp,cl_in,cln,clb,self.wsp.ncls*self.wsp.bin.n_bands)
        clout=np.reshape(cl1d,[self.wsp.ncls,self.wsp.bin.n_bands])

        return clout

def deprojection_bias(f1,f2,cls_guess) :
    if(len(cls_guess)!=f1.fl.nmaps*f2.fl.nmaps) :
        raise KeyError("Proposal Cell doesn't match number of maps")
    if(len(cls_guess[0])!=f1.fl.lmax+1) :
        raise KeyError("Proposal Cell doesn't match map resolution")
    cl1d=lib.comp_deproj_bias(f1.fl,f2.fl,cls_guess,len(cls_guess)*len(cls_guess[0]))
    cl2d=np.reshape(cl1d,[len(cls_guess),len(cls_guess[0])])

    return cl2d

def compute_coupled_cell(f1,f2) :
    if(f1.fl.nside!=f2.fl.nside) :
        raise KeyError("Fields must have same resolution")
    
    cl1d=lib.comp_pspec_coupled(f1.fl,f2.fl,f1.fl.nmaps*f2.fl.nmaps*(f1.fl.lmax+1))
    clout=np.reshape(cl1d,[f1.fl.nmaps*f2.fl.nmaps,f1.fl.lmax+1])

    return clout
    
def compute_decoupled_cell(f1,f2,b,cl_noise=None,cl_guess=None,workspace=None) :
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

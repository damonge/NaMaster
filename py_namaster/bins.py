import nmtlib as lib
import numpy as np

class NmtBin(object) :

    def __init__(self,nside,fname=None,nlb=None) :
        if((nlb==None) and (fname==None)) :
            raise KeyError("Must supply filename or bandpower width")

        if(nlb==None) :
            self.bin=lib.bins_read(fname,3*nside-1)
        if(fname==None) :
            self.bin=lib.bins_create(nlb,3*nside-1)
        self.lmax=3*nside-1

    def __del__(self) :
        lib.bins_free(self.bin)

    def get_n_bands(self) :
        return self.bin.n_bands

    def get_nell_list(self) :
        return lib.get_nell_list(self.bin,self.bin.n_bands)

    def get_ell_list(self,ibin) :
        return lib.get_ell_list(self.bin,ibin,lib.get_nell(self.bin,ibin))

    def get_weight_list(self,ibin) :
        return lib.get_weight_list(self.bin,ibin,lib.get_nell(self.bin,ibin))

    def get_effective_ells(self) :
        return lib.get_ell_eff(self.bin,self.bin.n_bands)

    def bin_cell(self,cls_in) :
        if(len(cls_in[0])!=self.lmax+1) :
            raise KeyError("Input Cl has wrong size")
        cl1d=lib.bin_cl(self.bin,cls_in,len(cls_in)*self.bin.n_bands)
        clout=np.reshape(cl1d,[len(cls_in),self.bin.n_bands])
        return clout

    def unbin_cell(self,cls_in) :
        if(len(cls_in[0])!=self.bin.n_bands) :
            raise KeyError("Input Cl has wrong size")
        cl1d=lib.unbin_cl(self.bin,cls_in,len(cls_in)*(self.lmax+1))
        clout=np.reshape(cl1d,[len(cls_in),self.lmax+1])
        return clout

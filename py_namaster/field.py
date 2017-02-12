import numpy as np
import nmtlib as lib

class NmtField(object) :

    def __init__(self,mask,maps,templates=None) :
        if((len(maps)!=1) and (len(maps)!=2)) :
            raise KeyError("Must supply 1 or 2 maps per field")
        if(templates!=None) :
            if((len(templates[0])!=1) and (len(templates[0])!=2)) :
                raise KeyError("Must supply 1 or 2 maps per field")

        if(templates==None) :
            self.fl=lib.field_alloc_new_notemp(mask,maps)
        else :
            self.fl=lib.field_alloc_new(mask,maps,templates)

    def __del__(self) :
        lib.field_free(self.fl)

    def get_maps(self) :
        maps=np.zeros([self.fl.nmaps,self.fl.npix])
        for imap in np.arange(self.fl.nmaps) :
            maps[imap,:]=lib.get_map(self.fl,imap,int(self.fl.npix))
        return maps

    def get_templates(self) :
        temp=np.zeros([self.fl.ntemp,self.fl.nmaps,self.fl.npix])
        for itemp in np.arange(self.fl.ntemp) :
            for imap in np.arange(self.fl.nmaps) :
                temp[itemp,imap,:]=lib.get_temp(self.fl,itemp,imap,int(self.fl.npix))
        return temp

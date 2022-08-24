#!/usr/bin/env python3
__all__ = ['read_dms']
import numpy as np

def read_dms(dmskey)
    tmp=np.fromfile(dmskey,dtype='<d')
    return tmp

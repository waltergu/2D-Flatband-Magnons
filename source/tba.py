import numpy as np
import HamiltonianPy.FreeSystem as TBA
from HamiltonianPy import *
from config import *

__all__=['tbaconstruct']

def tbaconstruct(name,parameters,lattice,terms,**karg):
    config=IDFConfig(priority=DEFAULT_FERMIONIC_PRIORITY,pids=lattice.pids,map=idfmap)
    tba=TBA.TBA(
        dlog=       'log',
        din=        'data',
        dout=       'result/tba',
        name=       '%s_%s'%(name,lattice.name),
        parameters= parameters,
        map=        parametermap,
        lattice=    lattice,
        config=     config,
        terms=      [term(**parameters) for term in terms],
        dtype=      np.complex128
        )
    return tba

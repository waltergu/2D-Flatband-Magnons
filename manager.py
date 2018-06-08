import mkl
import numpy as np
from HamiltonianPy import *
from source import *
from collections import OrderedDict
from fractions import Fraction

def tbatasks(name,parameters,lattice,terms,jobs=()):
    import HamiltonianPy.FreeSystem as TBA
    tba=tbaconstruct(name,parameters,lattice,terms)
    if 'EB' in jobs:
        if len(lattice.vectors)==2:
            tba.register(EB(name='EB',path=square_gxm(reciprocals=lattice.reciprocals,nk=100),run=TBA.TBAEB))
        elif len(lattice.vectors)==1:
            tba.register(EB(name='EB',path=KSpace(reciprocals=lattice.reciprocals,segments=[(-0.5,0.5)],end=True,nk=401),run=TBA.TBAEB))
        else:
            tba.register(EB(name='EB',run=TBA.TBAEB))
        tba.summary()

def fbfmtasks(name,parameters,lattice,terms,interactions,nk=50,scalefree=1.0,scaleint=1.0,jobs=()):
    import HamiltonianPy.FBFM as FB
    assert  len(lattice.vectors)==2
    ns,ne=len(lattice),len(lattice)/2
    basis=FB.FBFMBasis(BZ=FBZ(lattice.reciprocals,nks=(nk,nk)),filling=Fraction(ne,ns*2))
    fbfm=fbfmconstruct(name,parameters,basis,lattice,terms,interactions)
    if 'EB' in jobs:
        fbfm.register(FB.EB(name='EB%s'%nk,path='S:G-Y,Y-M,M-G',ne=nk**2,scalefree=scalefree,scaleint=scaleint,plot=True,run=FB.FBFMEB))
        fbfm.summary()
    if 'KPOS' in jobs:
        path=[fbfm.basis.BZ[pos] for pos in fbfm.basis.BZ.path(KMap(fbfm.lattice.reciprocals,'S:G-X,X-M,M-G'),mode='I')]
        pos=len(path)/6*5
        mode='D'
        fbfm.register(KPOS(name='KPOS%s'%nk,k=path[pos],ns=(0,1,3,4),plot=True,scalefree=0.001,mode=mode,run=FBFMKPOS))
        fbfm.summary()

if __name__=='__main__':
    mkl.set_num_threads(1)
    Engine.DEBUG=True
    Engine.MKDIR=False

    from ConfigParser import SafeConfigParser
    config=SafeConfigParser()
    config.read('input.cfg')
    jobs=config.get('controllers','jobs')

    # parameters
    parameters=OrderedDict()
    parameters['t1']=-1.0
    parameters['t2']=float(config.get('parameters','t2'))

    # tba
    m=int(config.get('controllers','m'))
    if 'CITBA0' in jobs: tbatasks(name1,parameters,S2xxy('1P-1P',nneighbour),[t1c,t2],jobs=['EB'])
    if 'CITBA1' in jobs: tbatasks(name1,parameters,S2xxy('%sO-1P'%m,nneighbour),[t1c,t2],jobs=['EB'])
    if 'TITBA0' in jobs: tbatasks(name2,parameters,S2xxy('1P-1P',nneighbour),[t1t,t2],jobs=['EB'])
    if 'TITBA1' in jobs: tbatasks(name2,parameters,S2xxy('%sO-1P'%m,nneighbour),[t1t,t2],jobs=['EB'])

    # fbfm
    nk=int(config.get('controllers','nk'))
    parameters['U']=float(config.get('parameters','U'))
    if 'CIEB0' in jobs: fbfmtasks(name1,parameters,S2xxy('1P-1P',nneighbour),[t1c,t2],[U],nk=nk,scalefree=1.0,scaleint=1.0,jobs=['EB'])
    if 'CIEB1' in jobs: fbfmtasks(name1,parameters,S2xxy('1P-1P',nneighbour),[t1c,t2],[U],nk=nk,scalefree=0.0,scaleint=1.0,jobs=['EB'])
    if 'CIEB2' in jobs: fbfmtasks(name1,parameters,S2xxy('1P-1P',nneighbour),[t1c,t2],[U],nk=nk,scalefree=1.0,scaleint=0.0,jobs=['EB'])
    if 'CIKPOS' in jobs: fbfmtasks(name1,parameters,S2xxy('1P-1P',nneighbour),[t1c,t2],[U],nk=nk,scalefree=1.0,scaleint=1.0,jobs=['KPOS'])

    if 'TIEB0' in jobs: fbfmtasks(name2,parameters,S2xxy('1P-1P',nneighbour),[t1t,t2],[U],nk=nk,scalefree=1.0,scaleint=1.0,jobs=['EB'])
    if 'TIEB1' in jobs: fbfmtasks(name2,parameters,S2xxy('1P-1P',nneighbour),[t1t,t2],[U],nk=nk,scalefree=0.0,scaleint=1.0,jobs=['EB'])
    if 'TIEB2' in jobs: fbfmtasks(name2,parameters,S2xxy('1P-1P',nneighbour),[t1t,t2],[U],nk=nk,scalefree=1.0,scaleint=0.0,jobs=['EB'])
    if 'TIKPOS' in jobs: fbfmtasks(name2,parameters,S2xxy('1P-1P',nneighbour),[t1t,t2],[U],nk=nk,scalefree=1.0,scaleint=1.0,jobs=['KPOS'])

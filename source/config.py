from HamiltonianPy import *
import numpy as np

__all__=['name1','name2','nneighbour','parametermap','idfmap','t1c','t1t','t2','U','S2xxy','KPOS','FBFMKPOS']

# The configs of the model
name1='2DCI'
name2='2DTI'
nneighbour=2

# parametermap
parametermap=None

# idfmap
idfmap=lambda pid: Fock(atom=0,norbital=1,nspin=2,nnambu=1)

# nearest neighbour pi-flux hopping indexpacks for Chern insulator
def cipi(bond):
    theta=azimuthd(bond.rcoord)
    s1,s2,phase=bond.spoint.pid.site%2,bond.epoint.pid.site%2,np.exp(1.0j*np.pi/4)
    if np.abs(theta)<RZERO or np.abs(theta-180)<RZERO:
        coeff=phase if (s1,s2)==(0,1) else (np.conjugate(phase) if (s1,s2)==(1,0) else 0.0)
    elif np.abs(theta-90)<RZERO or np.abs(theta-270)<RZERO:
        coeff=np.conjugate(phase) if (s1,s2)==(0,1) else (phase if (s1,s2)==(1,0) else 0.0)
    else:
        coeff=0.0
    return sigma0('sp')*coeff

# nearest neighbour pi-flux hopping indexpacks for topological insulator
def tipi(bond):
    theta=azimuthd(bond.rcoord)
    s1,s2,phase=bond.spoint.pid.site%2,bond.epoint.pid.site%2,np.exp(1.0j*np.pi/4)
    if np.abs(theta)<RZERO or np.abs(theta-180)<RZERO:
        coeff=phase if (s1,s2)==(0,1) else (np.conjugate(phase) if (s1,s2)==(1,0) else 0.0)
    elif np.abs(theta-90)<RZERO or np.abs(theta-270)<RZERO:
        coeff=np.conjugate(phase) if (s1,s2)==(0,1) else (phase if (s1,s2)==(1,0) else 0.0)
    else:
        coeff=0.0
    return FockPack(coeff,spins=(1,1))+FockPack(np.conjugate(coeff),spins=(0,0))

# next to nearset neighbour hopping amplitude
def nnnamplitude(bond):
    theta=azimuthd(bond.rcoord)
    s1,s2=bond.spoint.pid.site%2,bond.epoint.pid.site%2
    if np.abs(theta-45)<RZERO or np.abs(theta-225)<RZERO:
        result=-1.0 if (s1,s2)==(0,0) else (1.0 if (s1,s2)==(1,1) else 0.0)
    elif np.abs(theta-135)<RZERO or np.abs(theta-315)<RZERO:
        result=1.0 if (s1,s2)==(0,0) else (-1.0 if (s1,s2)==(1,1) else 0.0)
    else:
        result=0.0
    return result

# terms
t1c=lambda **parameters: Hopping('t1',parameters['t1'],neighbour=1,indexpacks=cipi)
t1t=lambda **parameters: Hopping('t1',parameters['t1'],neighbour=1,indexpacks=tipi)
t2=lambda **parameters: Hopping('t2',parameters['t2'],neighbour=2,amplitude=nnnamplitude)
U=lambda **parameters: Hubbard('U',parameters['U'])

# clusters
S2xxy=Square('S2xxy')

# App and method
class PPOS(POS):
    def __init__(self,k,path,ns=(0,1,2),scalefree=0.001,**karg):
        self.k=k
        self.path=path
        self.ns=np.array(ns)
        self.scalefree=scalefree

def FBFMPPOS(engine,app):
    import matplotlib.pyplot as plt
    import pickle as pk
    from scipy.linalg import eigh
    path=np.array(engine.basis.BZ.path(KMap(engine.lattice.reciprocals,app.path),mode='I'))
    matrix=engine.matrix(k=app.k,scalefree=app.scalefree)
    vs0=(np.abs(eigh(matrix,eigvals_only=False)[1][:,app.ns])**2)[path]
    matrix=engine.matrix(k=app.k,scalefree=1.0)
    vs1=(np.abs(eigh(matrix,eigvals_only=False)[1][:,app.ns])**2)[path]
    result=vs1-vs0
    name='%s_%s_%s'%(engine.tostr(),app.name,repr(app.k))
    if app.savedata:
        with open('%s/%s.pkb'%(engine.dout,name),'wb') as fout: pk.dump(result,fout,2)
    if app.plot:
        plt.title(name)
        plt.plot(range(len(path)),result)
        if app.show and app.suspend: plt.show()
        if app.show and not app.suspend: plt.pause(App.SUSPEND_TIME)
        if app.savefig: plt.savefig('%s/%s.png'%(engine.dout,name))
        plt.close()
    if app.returndata: return result

class KPOS(POS):
    def __init__(self,k,ns=(0,1,2),mode='FND',scalefree=0.001,**karg):
        self.k=k
        self.ns=np.array(ns)
        self.mode=mode
        self.scalefree=scalefree

def FBFMKPOS(engine,app):
    import matplotlib.pyplot as plt
    import itertools as it
    import pickle as pk
    from scipy.linalg import eigh
    shape=tuple(it.chain([len(app.ns)],engine.basis.BZ.contents.max(axis=0)+1))
    permutation=np.argsort((engine.basis.BZ+engine.basis.BZ.type((shape[1]/2,shape[2]/2))).sorted(history=True)[1])
    permutation=(engine.basis.BZ+engine.basis.BZ.type((shape[1]/2,shape[2]/2))).sorted(history=True)[1]
    if 'F' in app.mode or 'D' in app.mode:
        matrix=engine.matrix(k=app.k,scalefree=app.scalefree)
        vs0=(np.abs(eigh(matrix,eigvals_only=False)[1][:,app.ns])**2).T[:,permutation].reshape(shape).transpose([0,2,1])
    if 'N' in app.mode or 'D' in app.mode:
        matrix=engine.matrix(k=app.k,scalefree=1.0)
        vs1=(np.abs(eigh(matrix,eigvals_only=False)[1][:,app.ns])**2).T[:,permutation].reshape(shape).transpose([0,2,1])
    result=[vs0,vs1] if len(app.mode)>1 else vs0 if app.mode=='F' else vs1 if app.mode=='N' else vs1-vs0
    name='%s_%s_%s'%(engine.tostr(),app.name,repr(app.k))
    if app.savedata:
        with open('%s/%s.pkb'%(engine.dout,name),'wb') as fout: pk.dump(result,fout,2)
    if app.plot:
        plt.axis('equal')
        for i,ne in enumerate(app.ns):
           for tag in app.mode:
               plt.title('%s_%s%s'%(name,ne,tag))
               plt.colorbar(plt.pcolormesh(vs0[i,:,:] if tag=='F' else vs1[i,:,:] if tag=='N' else vs1[i,:,:]-vs0[i,:,:]))
               if app.show and app.suspend: plt.show()
               if app.show and not app.suspend: plt.pause(App.SUSPEND_TIME)
               if app.savefig: plt.savefig('%s/%s_%s%s.png'%(engine.dout,name,ne,tag))
               plt.close()
    if app.returndata: return result

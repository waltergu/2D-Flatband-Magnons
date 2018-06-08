import HamiltonianPy as HP
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.interpolate as si
import sys
import pdb

def lattice():
    plt.ion()
    fig=plt.figure()
    gs=plt.GridSpec(1,2)
    fig.subplots_adjust(left=0.02,right=0.98,top=0.99,bottom=0.0,hspace=0.0,wspace=0.08)

    ax=fig.add_subplot(gs[0,0])
    ax.axis('off')
    ax.axis('equal')
    ax.set_xlim(-0.1,2.1)
    ax.set_ylim(-0.1,2.1)
    P,t1,t2=np.array([0.0,0.0]),np.array([1.0,0.0]),np.array([0.0,1.0])
    lattice=HP.Lattice(name='2DFB',rcoords=HP.tiling(cluster=[P],vectors=[t1,t2],translations=it.product(xrange(3),xrange(3))),neighbours=2)
    for bond in lattice.bonds:
        p1,p2=bond.spoint,bond.epoint
        if bond.neighbour==0:
            ax.scatter(p1.rcoord[0],p1.rcoord[1],s=100,edgecolors='none',color='blue' if p1.pid.site%2==0 else 'red',zorder=3)
        elif bond.neighbour==1:
            ax.plot([p1.rcoord[0],p2.rcoord[0]],[p1.rcoord[1],p2.rcoord[1]],lw=2,color='black',zorder=2)
            theta,s1,s2=HP.azimuthd(bond.rcoord),p1.pid.site%2,p2.pid.site%2
            if np.abs(theta)<HP.RZERO or np.abs(theta-180)<HP.RZERO:
                coeff=1 if (s1,s2)==(0,1) else -1
            else:
                coeff=-1 if (s1,s2)==(0,1) else 1
            center,direction=(p1.rcoord+p2.rcoord)/2,bond.rcoord/10*coeff
            ax.annotate(s='',xy=center-direction/2,xytext=center+direction/2,arrowprops={'color':'black','linewidth':2,'arrowstyle':'->','zorder':3})
        else:
            theta,s1,s2=HP.azimuthd(bond.rcoord),p1.pid.site%2,p2.pid.site%2
            if np.abs(theta-45)<HP.RZERO or np.abs(theta-225)<HP.RZERO:
                coeff=-1 if (s1,s2)==(0,0) else 1
            else:
                coeff=1 if (s1,s2)==(0,0) else -1
            ax.plot([p1.rcoord[0],p2.rcoord[0]],[p1.rcoord[1],p2.rcoord[1]],lw=2,zorder=2,ls='--',color='green' if coeff==1 else 'purple')
    ax.fill_between(np.array([0.0,1.0]),np.array([1.0,2.0]),np.array([1.0,0.0]),alpha=0.7,zorder=0,color='grey')
    ax.fill_between(np.array([1.0,2.0]),np.array([2.0,1.0]),np.array([0.0,1.0]),alpha=0.7,zorder=0,color='grey')
    ax.text(-0.1,2.1,"(a)",ha='left',fontsize=18,color='black')

    ax=fig.add_subplot(gs[0,1])
    ax.axis('off')
    ax.axis('equal')
    ax.set_xlim(-1.5,1.8)
    ax.set_ylim(-1.5,1.8)
    ax.annotate(s='',xy=[-1.2,0],xytext=[1.5,0],arrowprops={'color':'black','linewidth':2,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[0,-1.2],xytext=[0,1.5],arrowprops={'color':'black','linewidth':2,'arrowstyle':'<-','zorder':3})
    p0,p1,p2,p3,p4=np.array([0.0,0.0]),np.array([-1.0,-1.0]),np.array([-1.0,1.0]),np.array([1.0,1.0]),np.array([1.0,-1.0])
    for sp,ep in [(p1,p2),(p2,p3),(p3,p4),(p4,p1),(p0,p3)]:
        ax.plot([sp[0],ep[0]],[sp[1],ep[1]],lw=2,color='black',zorder=1)
    ax.text(-0.1,0.1,'$\Gamma$',fontsize=18,va='bottom',ha='right',color='black')
    ax.text(0.9,-0.1,'$X_1$',fontsize=18,va='top',ha='right',color='black')
    ax.text(-0.9,-0.1,'$X_2$',fontsize=18,va='top',ha='left',color='black')
    ax.text(-0.1,1.05,'$Y_1$',fontsize=18,va='bottom',ha='right',color='black')
    ax.text(-0.1,-1.05,'$Y_2$',fontsize=18,va='top',ha='right',color='black')
    ax.text(1.0,1.05,'$M_1$',fontsize=18,va='bottom',ha='center',color='black')
    ax.text(-1.0,1.05,'$M_2$',fontsize=18,va='bottom',ha='center',color='black')
    ax.text(-1.0,-1.05,'$M_3$',fontsize=18,va='top',ha='center',color='black')
    ax.text(1.0,-1.05,'$M_4$',fontsize=18,va='top',ha='center',color='black')
    ax.text(1.45,-0.1,'$k_x$',fontsize=12,va='top',ha='right',color='black')
    ax.text(0.1,1.45,'$k_y$',fontsize=12,va='top',ha='left',color='black')
    ax.text(-1.5,1.8,'(b)',ha='left',va='top',fontsize=18,color='black')

    pdb.set_trace()
    plt.savefig('lattice.pdf')
    plt.close()

def phase():
    plt.ion()
    fig,axes=plt.subplots(nrows=1,ncols=2)
    fig.subplots_adjust(left=0.10,right=0.875,top=0.97,bottom=0.185,hspace=0.2,wspace=0.310)

    for i,ax in enumerate(axes):
        xs=np.array([0.40,0.45,0.50,0.55,0.60,0.65,1/np.sqrt(2),0.75,0.80])
        if i==0:
            ys=np.array([2.53,2.16,1.89,1.83,1.79,1.81,2.00,2.20,2.50])
        else:
            ys=np.array([2.53,2.14,1.82,1.78,1.77,1.76,1.82,2.07,2.46])
        X=np.linspace(xs.min(),xs.max(),10*len(xs))
        tck=si.splrep(xs,ys,k=3)
        Y=si.splev(X,tck,der=0)
        ax.plot(X,Y,lw=2,color='black',zorder=1)
        ax.minorticks_on()
        ax.set_xlim(0.4,0.8)
        ax.set_ylim(1.6,2.6)
        ax.set_xticks(np.linspace(0.4,0.8,5))
        ax.set_xticklabels(['0.4','0.5','0.6','0.7','0.8'])
        ax.set_yticks(np.linspace(1.4,2.6,7))
        for tick in ax.get_xticklabels():
            tick.set_fontsize(18)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(18)
        ax.set_xlabel("$t_2/t_1$",fontdict={'fontsize':22})
        ax.set_ylabel("$U/t_1$",fontdict={'fontsize':20})
        ax.text(0.60,1.60,"NFM",va='center',ha='center',fontsize=22,color='black')
        ax.text(0.60,2.20,"FM",va='center',ha='center',fontsize=22,color='black')
        if i==0:
            ax.scatter([0.490],[2.3],s=50,marker='*',edgecolors='red',color='red',alpha=0.5,zorder=4)
            ax.scatter([0.430],[2.3],s=50,marker='*',edgecolors='blue',color='blue',alpha=0.5,zorder=4)
        else:
            ax.scatter([0.650],[2.0],s=50,marker='*',edgecolors='red',color='red',alpha=0.5,zorder=4)
            ax.scatter([0.741],[2.0],s=50,marker='*',edgecolors='blue',color='blue',alpha=0.5,zorder=4)

        t2s,Us=np.linspace(0.4,0.8,101),np.linspace(1.4,2.6,121)
        ks=HP.square_gxm(nk=200).mesh('k')
        d1k=lambda kx,ky,t1: -2*np.sqrt(2)*t1*np.cos(kx/2)*np.cos(ky/2)
        d2k=lambda kx,ky,t1: -2*np.sqrt(2)*t1*np.sin(kx/2)*np.sin(ky/2)
        d3k=lambda kx,ky,t2: -2*t2*(np.cos(kx)-np.cos(ky))
        ek=lambda kx,ky,t1,t2: np.sqrt(d1k(kx,ky,t1)**2+d2k(kx,ky,t1)**2+d3k(kx,ky,t2)**2)
        eks=lambda t1,t2: np.array([ek(kx,ky,t1,t2) for kx,ky in ks])
        EKS=[eks(1.0,t2) for t2 in t2s]
        widths=np.array([EK.max()-EK.min() for EK in EKS])
        flatness=1/widths
        Xs=np.tensordot(t2s,np.ones(len(Us)),axes=0)
        Ys=np.tensordot(np.ones(len(t2s)),Us,axes=0)
        Zs=np.tensordot(flatness,Us,axes=0)
        pcolor=ax.pcolormesh(Xs,Ys,Zs,alpha=1.0,cmap='coolwarm')
        if i==1:
            cbarax=fig.add_axes([0.91,0.195,0.04,0.70])
            colorbar=fig.colorbar(pcolor,cax=cbarax)
            colorbar.ax.set_title('$U/W$',fontdict={'fontsize':20})

    pdb.set_trace()
    plt.savefig('phase.pdf')
    plt.close()

def fmcispectrum():
    import matplotlib.gridspec as mg
    plt.ion()

    fig=plt.figure()
    gs=mg.GridSpec(4,4,width_ratios=[1,1,1,1],height_ratios=[1,1,1.75,1.75])
    fig.subplots_adjust(left=0.11,right=0.98,top=0.98,bottom=0.06,hspace=0.8,wspace=0.9)

    t,U,nk=-0.49,2.3,52
    posx,posm=26,52

    ax=fig.add_subplot(gs[0:2,0:2])
    data=np.load('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_EB%s.npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
    ax.plot(data[:,0],data[:,1:3]/U,color='green',lw=2,zorder=3)
    ax.plot(data[:,0],data[:,3:]/U,color='grey',lw=2,alpha=0.3,zorder=1)
    data=np.load('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_EB%s(0.0,1.0).npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
    ax.plot(data[:,0],data[:,1:3]/U,color='green',ls='--',lw=1.5,alpha=0.5,zorder=4)
    ax.plot(data[:,0],data[:,3:5]/U,color='green',ls='--',lw=1.5,alpha=0.2,zorder=3)
    data=np.load('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_EB%s(1.0,0.0).npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
    ax.plot(data[:,0],data[:,+1]/U+0.5,color='blue',lw=1,zorder=2,alpha=0.9)
    ax.plot(data[:,0],data[:,-1]/U+0.5,color='blue',lw=1,zorder=2,alpha=0.9)
    ax.axvline(x=posx,ls='--',color='grey',lw=1.5)
    ax.axvline(x=posm,ls='--',color='grey',lw=1.5)
    ax.minorticks_on()
    ax.set_xlim(0,len(data)-1)
    ax.set_ylim(0.0,1.0)
    ax.set_xticks([0,posx,posm,len(data)-1])
    ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'])
    ax.set_yticks(np.linspace(0.0,1.0,3))
    ax.set_yticklabels(['0','0.5','1.0'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel('q',fontdict={'fontsize':18})
    ax.set_ylabel('$E/U$',fontdict={'fontsize':18})
    ax.annotate(s='',xy=[posx,0.11],xytext=[posx,0.22],arrowprops={'color':'black','linewidth':2,'arrowstyle':'->','zorder':3})
    ax.annotate(s='',xy=[posx/2*5,0.13],xytext=[posx/2*5,0.23],arrowprops={'color':'black','linewidth':2,'arrowstyle':'->','zorder':3})
    ax.text(posx,0.22,'$q_1$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(posx/2*5,0.23,'$q_2$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(len(data)-2,0.85,'(a)',ha='right',va='bottom',fontsize=18,color='black')

    ax=fig.add_subplot(gs[0:2,2:])
    data=np.loadtxt('../result/tba/2DCI_S2xxy(1P-1P)_-1.0_%s_TBA_EB.dat'%HP.decimaltostr(t))
    emin=data[:,1].min()
    ax.plot(data[:,0],(data[:,1]-emin)/U,lw=1.5,color='blue',zorder=1,label='up')
    ax.plot(data[:,0],(data[:,1]-emin)/U+0.5,lw=1.5,color='blue',ls='--',zorder=1,label='down')
    ax.axvline(x=100,ls='--',color='grey',lw=1,alpha=0.9)
    ax.axvline(x=200,ls='--',color='grey',lw=1,alpha=0.9)
    ax.minorticks_on()
    ax.set_xlim(0,300)
    ax.set_ylim(0.0,1.0)
    ax.set_xticks([0,100,200,300])
    ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'])
    ax.set_yticks([0.0,0.5,1.0])
    ax.set_yticklabels(['0.0','0.5','1.0'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel('k',fontdict={'fontsize':18})
    ax.set_ylabel('$(E-E_{min})/U$',fontdict={'fontsize':18})
    ax.text(299,0.85,'(c)',ha='right',va='bottom',fontsize=18,color='black')
    leg=ax.legend(loc='lower left',fancybox=True,shadow=False,prop={'size': 14})
    leg.get_frame().set_alpha(0.5)

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax=inset_axes(ax,width="70%",height=1.8,loc='upper left')
    ax.axis('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0.8)
    ax.set_xlim(-1.45,1.45)
    ax.set_ylim(-1.45,1.45)
    p1,p2,p3,p4=[-1.0,-1.0],[-1.0,1.0],[1.0,1.0],[1.0,-1.0]
    for (x1,y1),(x2,y2) in [(p1,p2),(p2,p3),(p3,p4),(p4,p1)]:
        ax.plot([x1,x2],[y1,y2],lw=2,color='black')
    ax.text(-0.04,+0.00,'$\Gamma$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+1.04,+1.00,'$M_1$',ha='left',va='center',fontsize=14,color='black')
    ax.text(-1.04,+1.00,'$M_2$',ha='right',va='center',fontsize=14,color='black')
    ax.text(-1.04,-1.00,'$M_3$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+1.04,-1.00,'$M_4$',ha='left',va='center',fontsize=14,color='black')
    ax.text(+1.04,+0.00,'$X_1$',ha='left',va='center',fontsize=14,color='black')
    ax.text(-1.04,+0.00,'$X_2$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+0.00,+0.96,'$Y_1$',ha='center',va='top',fontsize=14,color='black')
    ax.text(+0.00,-0.96,'$Y_2$',ha='center',va='bottom',fontsize=14,color='black')
    for x,y in ([0.0,0.0],[1.0,1.0],[-1.0,1.0],[-1.0,-1.0],[1.0,-1.0]):
        ax.scatter(x,y,s=50,edgecolors='none',color='purple',alpha=0.8,marker='o',zorder=3)
    for x,y in ([1.0,0.0],[-1.0,0.0],[0.0,1.0],[0.0,-1.0],[0.5,0.5],[-0.5,0.5],[-0.5,-0.5],[0.5,-0.5]):
        ax.scatter(x,y,s=50,edgecolors='none',color='red',alpha=0.8,marker='o',zorder=3)
    ax.text(+0.5,+0.54,'$O_1$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(-0.5,+0.54,'$O_2$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(-0.5,-0.54,'$O_3$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(+0.5,-0.54,'$O_4$',ha='center',va='bottom',fontsize=14,color='black')
    ax.annotate(s='',xy=[+0.0,+0.0],xytext=[+1.0,+0.0],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[-1.0,-1.0],xytext=[+0.0,-1.0],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.text(+0.5,+0.0,'$q_1$',ha='center',va='center',fontsize=14,color='black')
    ax.text(-0.5,-1.0,'$q_1$',ha='center',va='center',fontsize=14,color='black')
    ax.annotate(s='',xy=[+0.0,+0.0],xytext=[+0.5,+0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[-1.0,-1.0],xytext=[-0.5,-0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.text(+0.25,+0.25,'$q_2$',ha='center',va='center',fontsize=14,color='black')
    ax.text(-0.75,-0.75,'$q_2$',ha='center',va='center',fontsize=14,color='black')

    nk=40
    import pickle as pk
    import matplotlib.colorbar as mc
    igs=mg.GridSpecFromSubplotSpec(2,4,gs[2:,:],wspace=0.3,hspace=0.2)
    xs,ys=np.tile(np.array(xrange(nk+1)),nk+1),np.repeat(np.array(xrange(nk+1)),nk+1)
    with open('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_KPOS%s_kp(0,0).bin'%(HP.decimaltostr(t),HP.decimaltostr(U),nk),'rb') as fin:
        data=pk.load(fin)
    for i in xrange(4):
        result=np.zeros((nk+1,nk+1))
        result[:nk,:nk]=data[i,:,:]
        result[nk,:nk]=result[0,:nk]
        result[:nk,nk]=result[:nk,0]
        result[nk,nk]=result[0,0]
        result=result.reshape(-1)
        ax=fig.add_subplot(igs[i/2,i%2])
        ticks=[0,0.003] if i in (0,1) else [0,0.8]
        labels=['0','0.003'] if i in (0,1) else ['0','0.8']
        if i in (2,3):
           masks=np.abs(result)<(ticks[-1])*0.7
           permutation=np.concatenate([np.nonzero(masks)[0],np.nonzero(~masks)[0]])
        Xs,Ys,Result=(xs,ys,result) if i in (0,1) else (xs[permutation],ys[permutation],result[permutation])
        pcolor=ax.scatter(Xs,Ys,c=Result,vmin=ticks[0],vmax=ticks[-1],s=4,cmap=plt.cm.get_cmap("Reds"),marker='s',clip_on=False,zorder=2)
        cax,options=mc.make_axes(ax,location='top',shrink=0.6)
        cbar=fig.colorbar(pcolor,cax=cax,ticks=ticks,**options)
        cbar.set_ticklabels(labels)
        ax.axis('equal')
        ax.set_xlim(0,nk)
        ax.set_ylim(0,nk)
        ax.set_aspect('equal',adjustable='box')
        ax.set_xticks([0,nk/2,nk])
        ax.set_xticklabels(['$-1$','$0$','$1$'])
        ax.set_yticks([0,nk/2,nk])
        ax.set_yticklabels(['$-1$','$0$','$1$'])
        for tick in ax.get_xticklabels():
            tick.set_fontsize(14)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(14)
        if i in (2,3): ax.set_xlabel('$k_x/\pi$',fontdict={'fontsize':16})
        if i in (0,2): ax.set_ylabel('$k_y/\pi$',fontdict={'fontsize':16})
        for stag in ('top','bottom','left','right'):
            ax.spines[stag].set_alpha(0.3)
            ax.spines[stag].set_color('black')
        ax.tick_params(color='grey',grid_alpha=0.3)
        ax.text(-nk/3,nk+nk/3,'($b_%s$)'%(i+1),ha='left',va='top',fontsize=14,color='black')

    with open('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_KPOS%s_kp(%s,0).bin'%(HP.decimaltostr(t),HP.decimaltostr(U),nk,nk/2),'rb') as fin:
        data=pk.load(fin)
    for i in xrange(4):
        result=np.zeros((nk+1,nk+1))
        result[:nk,:nk]=data[i,:,:]
        result[nk,:nk]=result[0,:nk]
        result[:nk,nk]=result[:nk,0]
        result[nk,nk]=result[0,0]
        result=result.reshape(-1)
        ax=fig.add_subplot(igs[i/2,i%2+2])
        ticks=[-0.06,0.06] if i in (0,1) else [-0.06,0.06]
        labels=['-0.06','0.06'] if i in (0,1) else ['-0.06','0.06']
        if i in (2,3):
           masks=np.abs(result)<(ticks[-1])*0.7
           permutation=np.concatenate([np.nonzero(masks)[0],np.nonzero(~masks)[0]])
        Xs,Ys,Result=(xs,ys,result) if i in (0,1) else (xs[permutation],ys[permutation],result[permutation])
        pcolor=ax.scatter(Xs,Ys,c=Result,vmin=ticks[0],vmax=ticks[-1],s=4,cmap=plt.cm.get_cmap("coolwarm"),marker='s',clip_on=False,zorder=2)
        cax,options=mc.make_axes(ax,location='top',shrink=0.6)
        cbar=fig.colorbar(pcolor,cax=cax,ticks=ticks,**options)
        cbar.set_ticklabels(labels)
        ax.axis('equal')
        ax.set_xlim(0,nk)
        ax.set_ylim(0,nk)
        ax.set_aspect('equal',adjustable='box')
        ax.set_xticks([0,nk/2,nk])
        ax.set_xticklabels(['$-1$','$0$','$1$'])
        ax.set_yticks([0,nk/2,nk])
        ax.set_yticklabels(['$-1$','$0$','$1$'])
        for tick in ax.get_xticklabels():
            tick.set_fontsize(14)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(14)
        if i in (2,3): ax.set_xlabel('$k_x/\pi$',fontdict={'fontsize':16})
        for stag in ('top','bottom','left','right'):
            ax.spines[stag].set_alpha(0.3)
            ax.spines[stag].set_color('black')
        ax.tick_params(color='grey',grid_alpha=0.3)
        ax.text(-nk/3,nk+nk/3,'($d_%s$)'%(i+1),ha='left',va='top',fontsize=14,color='black')

    pdb.set_trace()
    plt.savefig('fmcispectrum.pdf')
    plt.close()

def pbcispectrum():
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    plt.ion()

    fig,axes=plt.subplots(nrows=1,ncols=2)
    fig.subplots_adjust(left=0.11,right=0.98,top=0.97,bottom=0.195,hspace=0.2,wspace=0.310)

    t,U,nk=-0.43,2.3,52
    posx,posm=26,52

    ax=axes[0]
    data=np.load('../result/fbfm/2DCI_S2xxy(1P-1P)_up_-1.0_%s_%s_FBFM_EB%s.npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
    ax.plot(data[:,0],data[:,1:3]/U,color='green',lw=2,zorder=3)
    ax.plot(data[:,0],data[:,3:]/U,color='grey',lw=3,alpha=0.3,zorder=1)
    ax.axvline(x=posx,ls='--',color='grey',lw=1.5)
    ax.axvline(x=posm,ls='--',color='grey',lw=1.5)
    ax.minorticks_on()
    ax.set_xlim(0,len(data)-1)
    ax.set_ylim(0.0,1.0)
    ax.set_xticks([0,posx,posm,len(data)-1])
    ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'])
    ax.set_yticks(np.linspace(0.0,1.0,3))
    ax.set_yticklabels(['0','0.5','1.0'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel('q',fontdict={'fontsize':18})
    ax.set_ylabel('$E/U$',fontdict={'fontsize':18})
    ax.annotate(s='',xy=[posx,0.01],xytext=[posx,0.12],arrowprops={'color':'black','linewidth':2,'arrowstyle':'->','zorder':3})
    ax.text(posx,0.12,'$q_1$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(len(data)-2,0.85,'(a)',ha='right',va='bottom',fontsize=18,color='black')

    ax=inset_axes(ax,width="30%",height=0.8,loc=2)
    ax.patch.set_alpha(0.85)
    ax.plot(data[0:posx+1,0],data[0:posx+1,1:3]/U,color='green',lw=2,zorder=3)
    ax.minorticks_on()
    ax.set_xlim(0,posx)
    ax.set_ylim(0,0.04)
    ax.yaxis.tick_right()
    ax.set_xticks([0,posx])
    ax.set_xticklabels(['$\Gamma$','$X$'])
    ax.set_yticks([0.0,0.02,0.04])
    ax.set_yticklabels(['0','0.02','0.04'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(16)
        tick.set_color('black')
    for tick in ax.get_yticklabels():
        tick.set_fontsize(14)
        tick.set_color('black')

    ax=axes[1]
    data=np.loadtxt('../result/tba/2DCI_S2xxy(1P-1P)_-1.0_%s_TBA_EB.dat'%HP.decimaltostr(t))
    emin=data[:,1].min()
    ax.plot(data[:,0],(data[:,1]-emin)/U,lw=1.5,color='blue',zorder=1,label='up')
    ax.plot(data[:,0],(data[:,1]-emin)/U+0.5,lw=1.5,color='blue',ls='--',zorder=1,label='down')
    ax.axvline(x=100,ls='--',color='grey',lw=1,alpha=0.9)
    ax.axvline(x=200,ls='--',color='grey',lw=1,alpha=0.9)
    ax.minorticks_on()
    ax.set_xlim(0,300)
    ax.set_ylim(0.0,1.0)
    ax.set_xticks([0,100,200,300])
    ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'])
    ax.set_yticks([0.0,0.5,1.0])
    ax.set_yticklabels(['0.0','0.5','1.0'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel('k',fontdict={'fontsize':18})
    ax.set_ylabel('$(E-E_{min})/U$',fontdict={'fontsize':18})
    ax.text(299,0.85,'(b)',ha='right',va='bottom',fontsize=18,color='black')
    leg=ax.legend(loc='lower left',fancybox=True,shadow=False,prop={'size': 14})
    leg.get_frame().set_alpha(0.5)

    ax=inset_axes(ax,width="70%",height=1.8,loc='upper left')
    ax.axis('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_alpha(0.8)
    ax.set_xlim(-1.45,1.45)
    ax.set_ylim(-1.45,1.45)
    p1,p2,p3,p4=[-1.0,-1.0],[-1.0,1.0],[1.0,1.0],[1.0,-1.0]
    for (x1,y1),(x2,y2) in [(p1,p2),(p2,p3),(p3,p4),(p4,p1)]:
        ax.plot([x1,x2],[y1,y2],lw=2,color='black')
    ax.text(-0.04,+0.00,'$\Gamma$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+1.04,+1.00,'$M_1$',ha='left',va='center',fontsize=14,color='black')
    ax.text(-1.04,+1.00,'$M_2$',ha='right',va='center',fontsize=14,color='black')
    ax.text(-1.04,-1.00,'$M_3$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+1.04,-1.00,'$M_4$',ha='left',va='center',fontsize=14,color='black')
    ax.text(+1.04,+0.00,'$X_1$',ha='left',va='center',fontsize=14,color='black')
    ax.text(-1.04,+0.00,'$X_2$',ha='right',va='center',fontsize=14,color='black')
    ax.text(+0.00,+0.96,'$Y_1$',ha='center',va='top',fontsize=14,color='black')
    ax.text(+0.00,-0.96,'$Y_2$',ha='center',va='bottom',fontsize=14,color='black')
    for x,y in ([0.0,0.0],[1.0,1.0],[-1.0,1.0],[-1.0,-1.0],[1.0,-1.0]):
        ax.scatter(x,y,s=50,edgecolors='none',color='purple',alpha=0.8,marker='o',zorder=3)
    for x,y in ([1.0,0.0],[-1.0,0.0],[0.0,1.0],[0.0,-1.0],[0.5,0.5],[-0.5,0.5],[-0.5,-0.5],[0.5,-0.5]):
        ax.scatter(x,y,s=50,edgecolors='none',color='red',alpha=0.8,marker='o',zorder=3)
    ax.text(+0.5,+0.54,'$O_1$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(-0.5,+0.54,'$O_2$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(-0.5,-0.54,'$O_3$',ha='center',va='bottom',fontsize=14,color='black')
    ax.text(+0.5,-0.54,'$O_4$',ha='center',va='bottom',fontsize=14,color='black')
    ax.annotate(s='',xy=[+0.0,+0.0],xytext=[+1.0,+0.0],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[-1.0,-1.0],xytext=[+0.0,-1.0],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':3})
    ax.text(+0.5,+0.0,'$q_1$',ha='center',va='center',fontsize=14,color='black')
    ax.text(-0.5,-1.0,'$q_1$',ha='center',va='center',fontsize=14,color='black')

    pdb.set_trace()
    plt.savefig('pbcispectrum.pdf')
    plt.close()

def tispectrum():
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    plt.ion()
    fig,axes=plt.subplots(nrows=2,ncols=2)
    fig.subplots_adjust(left=0.10,right=0.98,top=0.98,bottom=0.110,hspace=0.1,wspace=0.27)

    ts,U,nk=[0.65,0.741],2.0,52
    posx,posm=26,52
    for i,(t,tag) in enumerate(zip(ts,['(a)','(b)'])):
        ax=axes[i,0]
        data=np.load('../result/fbfm/2DTI_S2xxy(1P-1P)_up_-1.0_-%s_%s_FBFM_EB%s.npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
        ax.plot(data[:,0],data[:,1:3]/U,color='green',lw=2,zorder=3)
        ax.plot(data[:,0],data[:,3:]/U,color='grey',lw=3,alpha=0.3,zorder=1)
        data=np.load('../result/fbfm/2DTI_S2xxy(1P-1P)_up_-1.0_-%s_%s_FBFM_EB%s(0.0,1.0).npz'%(HP.decimaltostr(t),HP.decimaltostr(U),nk))['data']
        ax.plot(data[:,0],data[:,1:3]/U,color='green',ls='--',lw=1.5,alpha=0.5,zorder=4)
        ax.plot(data[:,0],data[:,3:5]/U,color='green',ls='--',lw=1.5,alpha=0.2,zorder=3)
        ax.axvline(x=posx,ls='--',color='grey',lw=1.5)
        ax.axvline(x=posm,ls='--',color='grey',lw=1.5)
        ax.minorticks_on()
        ax.set_xlim(0,len(data)-1)
        ax.set_ylim(0.0,0.6)
        ax.set_xticks([0,posx,posm,len(data)])
        ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'] if i==1 else ['']*4)
        ax.set_yticks(np.linspace(0.0,1.0,3))
        ax.set_yticklabels(['0','0.5','1.0'])
        for tick in ax.get_xticklabels():
            tick.set_fontsize(18)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(18)
        if i==1:ax.set_xlabel('q',fontdict={'fontsize':18})
        ax.set_ylabel('$E/U$',fontdict={'fontsize':18})
        ax.text(len(data)-2,0.85,tag,ha='right',va='bottom',fontsize=18,color='black')
        ystart=0.08 if i==0 else 0.01
        ax.annotate(s='',xy=[posx/2*5,ystart],xytext=[posx/2*5,ystart+0.11],arrowprops={'color':'black','linewidth':2,'arrowstyle':'->','zorder':3})
        ax.text(posx/2*5,ystart+0.11,'$q_1$',ha='center',va='bottom',fontsize=14,color='black')

    for i,(t,tag) in enumerate(zip(ts,['(c)','(d)'])):
        data=np.loadtxt('../result/tba/2DTI_S2xxy(1P-1P)_-1.0_-%s_TBA_EB.dat'%HP.decimaltostr(t))
        ax=axes[i,1]
        emin=data[:,1].min()
        ax.plot(data[:,0],(data[:,1]-emin)/U,lw=1.5,color='blue',zorder=1,label='up')
        ax.plot(data[:,0],(data[:,1]-emin)/U+0.5,lw=1.5,color='blue',ls='--',zorder=1,label='down')
        ax.axvline(x=100,ls='--',color='grey',lw=1,alpha=0.9)
        ax.axvline(x=200,ls='--',color='grey',lw=1,alpha=0.9)
        ax.minorticks_on()
        ax.set_xlim(0,300)
        ax.set_ylim(0.0,1.0)
        ax.set_xticks([0,100,200,300])
        ax.set_xticklabels(['$\Gamma$','$X_1$','$M_1$','$\Gamma$'] if i==1 else ['']*4)
        ax.set_yticks([0.0,0.5,1.0])
        ax.set_yticklabels(['0.0','0.5','1.0'])
        for tick in ax.get_xticklabels():
            tick.set_fontsize(18)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(18)
        if i==1: ax.set_xlabel('k',fontdict={'fontsize':18})
        ax.set_ylabel('$(E-E_{min})/U$',fontdict={'fontsize':18})
        ax.text(299,0.85,tag,ha='right',va='bottom',fontsize=18,color='black')
        leg=ax.legend(loc='lower left',fancybox=True,shadow=False,prop={'size': 14})
        leg.get_frame().set_alpha(0.8)

        ax=inset_axes(ax,width="70%",height=1.8,loc='upper left')
        ax.axis('equal')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.patch.set_alpha(0.8)
        ax.set_xlim(-1.45,1.45)
        ax.set_ylim(-1.45,1.45)
        p1,p2,p3,p4=[-1.0,-1.0],[-1.0,1.0],[1.0,1.0],[1.0,-1.0]
        for (x1,y1),(x2,y2) in [(p1,p2),(p2,p3),(p3,p4),(p4,p1)]:
            ax.plot([x1,x2],[y1,y2],lw=2,color='black')
        ax.text(-0.04,+0.00,'$\Gamma$',ha='right',va='center',fontsize=14,color='black')
        ax.text(+1.04,+1.00,'$M_1$',ha='left',va='center',fontsize=14,color='black')
        ax.text(-1.04,+1.00,'$M_2$',ha='right',va='center',fontsize=14,color='black')
        ax.text(-1.04,-1.00,'$M_3$',ha='right',va='center',fontsize=14,color='black')
        ax.text(+1.04,-1.00,'$M_4$',ha='left',va='center',fontsize=14,color='black')
        ax.text(+1.04,+0.00,'$X_1$',ha='left',va='center',fontsize=14,color='black')
        ax.text(-1.04,+0.00,'$X_2$',ha='right',va='center',fontsize=14,color='black')
        ax.text(+0.00,+0.96,'$Y_1$',ha='center',va='top',fontsize=14,color='black')
        ax.text(+0.00,-0.96,'$Y_2$',ha='center',va='bottom',fontsize=14,color='black')
        for x,y in ([0.0,0.0],[1.0,1.0],[-1.0,1.0],[-1.0,-1.0],[1.0,-1.0],[1.0,0.0],[-1.0,0.0],[0.0,1.0],[0.0,-1.0]):
            ax.scatter(x,y,s=50,edgecolors='none',color='purple',alpha=0.8,marker='o',zorder=3)
        pos1,pos2=data[0:100,1].argmax()*1.0/100,data[100:200,1].argmax()*1.0/100
        for x,y in ([pos1,0.0],[-pos1,0.0],[0.0,pos1],[0.0,-pos1],
                    [1.0,pos2],[1.0,-pos2],[-1.0,pos2],[-1.0,-pos2],
                    [pos2,1.0],[-pos2,1.0],[pos2,-1.0],[-pos2,-1.0],
                    [0.5,0.5],[-0.5,0.5],[-0.5,-0.5],[0.5,-0.5]):
            ax.scatter(x,y,s=50,edgecolors='none',color='red',alpha=0.8,marker='o',zorder=3)
        ax.text(+0.5,+0.54,'$O_1$',ha='center',va='bottom',fontsize=14,color='black')
        ax.text(-0.5,+0.54,'$O_2$',ha='center',va='bottom',fontsize=14,color='black')
        ax.text(-0.5,-0.54,'$O_3$',ha='center',va='bottom',fontsize=14,color='black')
        ax.text(+0.5,-0.54,'$O_4$',ha='center',va='bottom',fontsize=14,color='black')
        if i==0:
            ax.annotate(s='',xy=[+0.0,+0.0],xytext=[+0.5,+0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':2,'alpha':0.7})
            ax.annotate(s='',xy=[-1.0,-1.0],xytext=[-0.5,-0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':2,'alpha':0.7})
            ax.text(+0.25,+0.25,'$q_1$',ha='center',va='center',fontsize=14,color='black')
            ax.text(-0.75,-0.75,'$q_1$',ha='center',va='center',fontsize=14,color='black')
        else:
            ax.annotate(s='',xy=[-1.0,+0.0],xytext=[-0.5,+0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':2,'alpha':0.7})
            ax.annotate(s='',xy=[+0.0,-1.0],xytext=[+0.5,-0.5],arrowprops={'color':'green','linewidth':2.5,'arrowstyle':'<-','zorder':2,'alpha':0.7})
            ax.text(-0.75,+0.25,'$q_1$',ha='center',va='center',fontsize=14,color='black')
            ax.text(+0.25,-0.75,'$q_1$',ha='center',va='center',fontsize=14,color='black')

    pdb.set_trace()
    plt.savefig('tispectrum.pdf')
    plt.close()

if __name__=='__main__':
    for arg in sys.argv:
        if arg in ('1','all'): lattice()
        if arg in ('2','all'): phase()
        if arg in ('3','all'): fmcispectrum()
        if arg in ('4','all'): pbcispectrum()
        if arg in ('5','all'): tispectrum()

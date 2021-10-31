#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 12:20:19 2021

@author: luzhiwei
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 10:33:57 2021

@author: luzhiwei
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 19:39:00 2021

@author: luzhiwei
"""

import find_neighbor
import numpy as np
import math
import sys
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
global Nx,Ny,Nz,rc,pos_file
Nx=1
Ny=1
Nz=1
mi=1
J1=0.001
Niter_midpoint = 100
error_midpoint = 1E-8
Bi=np.array([0,0,-5])*4.24726E-05  
dt=0.1   #fs
step=30000
pos_file='POSCAR2'
r,latt_s=find_neighbor.expand_cell(pos_file,Nx,Ny,Nz)
Magmom_length=1
damping_onsite=0.2
non_scale=0.1
gamma=41.39
# damping_offsite=0.05
NN,NL1=find_neighbor.find_neighbor_pbc(pos_file,Nx,Ny,Nz)
# find_neighbor.plot_atom(r)

def gen_rand_spinlattice():
    spinLattice = np.zeros((len(r),3),np.float64)
    mySum = 1;
    randomVector = []
    randomSum = 0  
    for ii in range(len(r)): 
        mySum = 1;
        randomVector = []
        randomSum = 0  
        np.random.seed(ii+1)
        for i in range(3):           
            randomNumber = np.random.random()
            randomVector.append(randomNumber)
            randomSum  += randomNumber
        coef = mySum/randomSum
        myNewList = [j * coef for j in randomVector]
        spinLattice[ii,:]=np.array(myNewList)      
    spinLattice_final=np.sqrt(spinLattice)*Magmom_length
    return spinLattice_final

def gen_spec_spinlattice():
    spinLattice = np.zeros((len(r),3),np.float64)
    spinLattice[0,:]=np.array([0,-0.2,0.8])
    spinLattice[1,:]=np.array([0,0.2,0.8])
    # spinLattice[2,:]=np.array([1,0.01,0])
    # print(spinLattice)
    return spinLattice

def Hamitonian(spin):
    Hamitonian_Hei=0
    # Hamitonian_Zee=np.ones((len(r),3))
    # Hamitonian_Zee=-Bi*Hamitonian_Zee
    for ii in np.arange(len(r)):
        for jj in np.arange(NN[ii].astype('int32')):
            Hamitonian_Hei=Hamitonian_Hei+np.sum((-J1*spin[ii,:]*spin[NL1[ii,jj].astype('int32'),:]))
        Hamitonian_Hei=Hamitonian_Hei-np.sum(Bi*spin[ii,:])
    Hamitonian_Hei=1/len(r)*Hamitonian_Hei
    # Hamitonian=Hamitonian_Hei*2         
    return Hamitonian_Hei

def find_diag_nonlocal_damping():
    damping=np.zeros((3,3,2))
    pair2index_damping=np.ones((len(r),len(r)))*2
    for ii in np.arange(len(r)):
        for jj in np.arange(len(r)):
            if ii==jj:
                pair2index_damping[ii,jj]=0
        for kk in np.arange(int(NN[ii])):
            pair2index_damping[ii,int(NL1[ii,kk])]=1
    damping[0,0,0] = damping_onsite
    damping[1,1,0] = damping_onsite
    damping[2,2,0] = damping_onsite
    damping[0,0,1] = damping_onsite*non_scale
    damping[1,1,1] = damping_onsite*non_scale
    damping[2,2,1] = damping_onsite*non_scale
    return pair2index_damping,damping

def find_Beff(delta,s_new):
    B_eff=np.zeros((len(r),3))
    B_eff2 = Bi
    for ii in np.arange(len(r)):
        B_eff1=np.zeros((1,3))
        for jj in np.arange(int(NN[ii])):
            B_eff1=B_eff1+J1*s_new[int(NL1[ii,jj]),:]       #######少了*2
        B_eff[ii,:]=-gamma * (B_eff1 + B_eff2)
    pair2index_damping,damping=find_diag_nonlocal_damping()
    B_damp=np.zeros((len(r),3))
    for iatom in np.arange(len(r)):
        B_eff3=np.zeros((1,3))
        for jatom in np.arange(len(r)):
            jp=int(pair2index_damping[iatom,jatom])
            if jp!=2:
                mxB=np.cross(s_new[jatom,:],B_eff[jatom,:])
                dmdt=(delta[jatom,:]/dt)
                B_eff3=B_eff3+np.dot(damping[:,:,jp],dmdt)
        B_damp[iatom,:]=B_eff3    
    return B_eff,B_damp

       
def spin_dyn(t,s_new):
    delta=np.zeros((len(r),3))
    tmp_delta=np.zeros((len(r),3))
    s_mid=np.zeros((len(r),3))
    s_update=np.zeros((len(r),3))
    tmp_s=np.zeros((1,3))
    mxB=np.zeros((1,3))
    norm_s=0
    for ii in np.arange(Niter_midpoint):
        converge=True
        for iatom in np.arange(len(r)):
            tmp_s  = ((delta[iatom,:]+2*s_new[iatom,:])/2).reshape((1,3))
            norm_s = np.linalg.norm(tmp_s)
            s_mid[iatom,:] = mi/norm_s*tmp_s
        B_eff,B_damp=find_Beff(delta,s_mid)    
        for iatom in np.arange(len(r)):
            mxB=np.cross(s_mid[iatom,:],B_eff[iatom,:]+B_damp[iatom,:])
            delta[iatom,:]=dt*mxB
        for iatom in np.arange(len(r)):
            if not np.linalg.norm(tmp_delta[iatom,:]-delta[iatom,:])< error_midpoint:
                converge=False
                tmp_delta=delta
                break
        if converge:
            break
    if ii== Niter_midpoint:
        print('Warning')
    for iatom in np.arange(len(r)):
        tmp_s=delta[iatom,:]+s_new[iatom,:]
        norm_s=np.linalg.norm(tmp_s)
        s_update[iatom,:]=mi/norm_s*tmp_s
   
    return s_update


def spin_update(s_update):
    s_new=s_update
    return s_new

def plotMagnetization(mx, my, mz):
	fig, ax = plt.subplots(figsize=(10,10))
	x=np.linspace(0,dt*step/1000,step)
	plt.plot(x,mx, label = 'mx')
	plt.plot(x,my, label = 'my')
	plt.plot(x,mz, label = 'mz')
	plt.title('Magnetic Moment vs Time(dt=%f fs)  ' %dt)
	plt.xlabel('t [ps]',fontdict={'weight': 'normal', 'size': 13})
	plt.ylabel('M [$\mu_B$]',fontdict={'weight': 'normal', 'size': 13})
	plt.legend(loc=2)


def update_quiver(num,Q,s_xyz):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    U=s_xyz[:,0,num]
    V=s_xyz[:,1,num]

    Q.set_UVC(U,V)
    
    return Q,


def main():
    global s_xyz
    file = open('moment_tot.txt','w')
    file.write('#_time (fs)'+' '+ '          '+'E'+'              '+'|M|'+'              '+'Mx' +'              '+ 'My'+'              '+'Mz')
    # file2 = open('moment.txt','w')
    # file2.write( '#_time (fs), site, |m|,m%x, m%y, m%z')
    # file3 = open('beff.txt','w')
    # file3.write( '#_time (fs), site, Bexch%x, Bexch%y, Bexch%z, Bext%x, Bext%y, Bext%z')
    # file4 = open('bdamp.txt','w')
    # file4.write(' #_time (fs), site, Bdamp%x, Bdamp%y, Bdamp%z')
    time=0
    # s_new=gen_spec_spinlattice()
    s_new= gen_rand_spinlattice()
    mx=np.zeros((step,1));my=np.zeros((step,1));mz=np.zeros((step,1))
   
    s_xyz=np.zeros((len(r),3,step))
   
    for t in np.arange(step):
        s_avg=np.sum(s_new,axis=0)/len(r)
        s_xyz[:,:,t]=s_new
        mx[t,0]=s_avg[0]
        my[t,0]=s_avg[1]
        mz[t,0]=s_avg[2]  
        s_update=spin_dyn(t,s_new)
        M=1/len(r)*np.sum(s_update,axis=0)
        E=Hamitonian(s_update)
        s_new=spin_update(s_update)
        if t%2000==0 or t==0:
            file.write('\n'+str('%.5f' %time)+'         '+str('%.8f'% E)+'         '+str('%.8f'% np.linalg.norm(M))+'         '+str('%.8f'% M[0])+'         '+str('%.8f'% M[1])+'         '+str('%.8f'% M[2]))
        time=t*dt/1000
    return s_xyz
    
    # file2.close()
    # file3.close()
    # file4.close()
        

if __name__=="__main__":  
    s_xyz=main()
    fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
    X=r[:,0]
    Y=r[:,1]
    Z=r[:,2]
    U = s_xyz[:,0,0]
    V = s_xyz[:,1,0]
    W = s_xyz[:,2,0]
    plt.xlim(-2,5) 
    plt.ylim(-2,5)
    Q=ax.quiver(X,Y,Z,U,V,W,length=0.05,color='b', normalize=True,arrow_length_ratio=0.3,linewidths=2)
    
    
    def update_quiver(num):
        """updates the horizontal and vertical vector components by a
        fixed increment on each frame
        """
        global Q
        Q.remove()
        
        U = s_xyz[:,0,num]
        V = s_xyz[:,1,num]
        W = s_xyz[:,2,num]
        Q=ax.quiver(X,Y,Z,U,V,W,length=0.5,color='b', normalize=True,arrow_length_ratio=0.3,linewidths=2)
        return Q,
    
    
    ani = animation.FuncAnimation(fig, update_quiver, frames=np.arange(0,step,50),interval=1, blit=False)
    fig.tight_layout()
    f='/Users/luzhiwei/Research/damping_python/aa.gif'
    ani.save(f, writer='imagemagick', fps=10)
   






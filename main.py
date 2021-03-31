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
from matplotlib.animation import FuncAnimation
global Nx,Ny,Nz,rc,pos_file
Nx=1
Ny=1
Nz=1
J1=1
Bi=np.array([-1,0,0])  #Tesla
dt=0.0001  #fs
step=10000
pos_file='Dimer'
r,latt_s=find_neighbor.expand_cell(pos_file,Nx,Ny,Nz)
Magmom_length=1   #uB
damping_onsite=0.01  
# damping_offsite=0.05
NN,NL1=find_neighbor.find_neighbor_npbc(pos_file,Nx,Ny,Nz)
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
    # print(np.linalg.norm(spinLattice_final[1,:]))
    return spinLattice_final

def gen_spec_spinlattice():
    spinLattice = np.zeros((len(r),3),np.float64)
    spinLattice[0,:]=np.array([1,-0.01,0])
    spinLattice[1,:]=np.array([1,0.02,0])
    spinLattice[2,:]=np.array([1,0.01,0])
    # print(spinLattice)
    return spinLattice

def Hamitonian(spin):
    Hamitonian_Hei=np.zeros((len(r),3))
    # Hamitonian_Zee=np.ones((len(r),3))
    # Hamitonian_Zee=-Bi*Hamitonian_Zee
    for ii in np.arange(len(r)):
        for jj in np.arange(NN[ii].astype('int32')):
            Hamitonian_Hei[ii,:]=Hamitonian_Hei[ii,:]+(-J1*spin[ii,:]*spin[NL1[ii,jj].astype('int32'),:])
    Hamitonian=Hamitonian_Hei          
    return Hamitonian


def find_Beff(spin):
    H=Hamitonian(spin)
    B_eff=np.divide(H,spin,np.zeros_like(H,dtype=np.float64),where=spin!=0)-Bi
    damping= find_diag_nonlocal_damping()
    # print(damping)
    B_eff_damp=np.zeros((len(r),3))
    for ii in np.arange(len(r)):
        for jj in np.arange(len(r)):
            B_eff_damp[ii,:]= B_eff_damp[ii,:]+((np.dot((damping[ii,jj,:,:]),derivate(spin[jj,:],B_eff[jj,:]).reshape(3,1))).reshape(1,3)/np.linalg.norm((spin[jj,:])))

    B_eff_total=B_eff+B_eff_damp
    return B_eff_total



def find_onsite_damping():
    damping=np.zeros((len(r),len(r),3,3))
    for ii in np.arange(len(r)):
        damping[ii,ii,:,:]=np.identity(3)
    return damping*damping_onsite

def find_differ_onsite_damping():
    damping=np.zeros((len(r),len(r),3,3))
    for ii in np.arange(len(r)):
        damping[ii,ii,:,:]=np.identity(3)
    for ii in np.arange(len(r)):
        for jj in np.arange(3):
            damping[ii,ii,jj,jj]=(jj+1)
    return damping*damping_onsite

def find_full_onsite_damping():
    damping=np.zeros((len(r),len(r),3,3))
    for ii in np.arange(len(r)):
        damping[ii,ii,:,:]=np.array([[1,0.1,0.1],[0.1,1,0.1],[0.1,0.1,1]])
    return damping*damping_onsite

def find_full_nonlocal_damping():
    damping=np.zeros((len(r),len(r),3,3))
    for ii in np.arange(len(r)):
        for jj in np.arange(len(r)):
            if ii ==jj:
                damping[ii,jj,:,:]=np.identity(3)
            else:
                damping[ii,jj,:,:]=np.array([[1,0.1,0.1],[0.1,1,0.1],[0.1,0.1,1]])*0.1              
    return damping*damping_onsite

def find_diag_nonlocal_damping():
    damping=np.zeros((len(r),len(r),3,3))
    for ii in np.arange(len(r)):
        for jj in np.arange(len(r)):
            if ii ==jj:
                damping[ii,jj,:,:]=np.identity(3)
            else:
                damping[ii,jj,:,:]=np.identity(3)*0.1              
    return damping*damping_onsite
def giromagneticRatio():
    return -41.39

def derivate(S,H):
    return giromagneticRatio()*(np.cross(S.T,H.T))


def update_spin(spin, n):
    old_spin_length=np.linalg.norm(spin,axis=1)
    B_eff=find_Beff(spin)
    result = (derivate(spin[n,:], B_eff[n,:]))
    spin[n,0] = spin[n,0] + dt*result[0]
    spin[n,1] = spin[n,1] + dt*result[1]
    spin[n,2] = spin[n,2] + dt*result[2]
    new_spin_length=np.linalg.norm(spin,axis=1)
    spin_final=spin.T*(old_spin_length/new_spin_length)
    return spin_final.T

def plotMagnetization(mx, my, mz):
	fig, ax = plt.subplots(figsize=(10,10))
	x=np.linspace(0,dt*step,step)
	plt.plot(x,mx, label = 'mx')
	plt.plot(x,my, label = 'my')
	plt.plot(x,mz, label = 'mz')
	plt.title('Magnetic Moment vs Time(dt=%f fs)  ' %dt)
	plt.xlabel('t [fs]',fontdict={'weight': 'normal', 'size': 13})
	plt.ylabel('M [$\mu_B$]',fontdict={'weight': 'normal', 'size': 13})
	plt.legend(loc=2)
# 	plt.show()

def plotPositions():
    fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
    x,y,z= np.array((r[:,0],r[:,1], r[:,2]))
    u=spinLattice[:,0]
    v=spinLattice[:,1]
    w=spinLattice[:,2]
    plt.xlim(-2, 5) 
    plt.ylim(-2, 5)
    ax.set_zlim(0,3)
    ax.quiver(x, y,z, u, v, w, color='b')
    plt.show()
    
    
# mx = np.zeros((len(r),step))
# my = np.zeros((len(r),step))
# mz = np.zeros((len(r),step))
spinLattice=gen_rand_spinlattice()
# spinLattice=gen_spec_spinlattice()
# plotPositions()
# for t in range(step):
#  	for s in range(len(r)):
#  	 	spin = update_spin(spinLattice, s)
#  	 	mx[s][t] = spin[s,0]
#  	 	my[s][t] = spin[s,1]
#  	 	mz[s][t] = spin[s,2]

# Mx=mx.mean(axis=0)
# My=my.mean(axis=0)
# Mz=mz.mean(axis=0)
# plotMagnetization(Mx,My,Mz)
# plt.savefig('dt_%f.png')
# plotPositions()
# 

#######Animation
fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
x,y,z= np.array((r[:,0],r[:,1], r[:,2]))
u=spinLattice[:,0]
v=spinLattice[:,1]
w=spinLattice[:,2]
plt.xlim(-2, 5) 
plt.ylim(-2, 5)
Q=ax.quiver(x, y,z, u, v, w, color='b')
# plt.show()
##animation 3D
def rotate(t):
    for x in range(len(r)):
        spin=update_spin(spinLattice,x)
    return spin

def update(step):
    global Q
    spin=rotate(step)
    sx = spin[:,0]
    sy = spin[:,1]
    sz = spin[:,2]   
    Q.remove()
    Q=ax.quiver(x, y, z, sx,sy,sz, length=0.8,color='b', normalize=True)
    
    # Q=ax.quiver(x, y, z, sx,sy,sz, length=0.8,color='b', normalize=True)
    # Q.set_UVC(sx,sy,sz)
    return Q,

ani = FuncAnimation(fig, update,blit=False, frames=step,interval=0.000000000000001,cache_frame_data=False)
# fig.tight_layout()
plt.show()
# ani.save('move1.gif', writer='imagemagick', fps=10)






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
Bi=np.array([10,0,0])
dt=0.0005
step=10000
pos_file='POSCAR_Fe'
r,latt_s=find_neighbor.expand_cell(pos_file,Nx,Ny,Nz)
Magmom_length=2
damping_onsite=0.5
damping_offsite=0.05
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
        for i in range(3):
            randomNumber = np.random.random()
            randomVector.append(randomNumber)
            randomSum  += randomNumber
        coef = mySum/randomSum
        myNewList = [j * coef for j in randomVector]
        spinLattice[ii,:]=np.array(myNewList)
    spinLattice_final=np.sqrt(spinLattice)*Magmom_length
    return spinLattice_final

def Hamitonian(spin):
    
    Hamitonian_Hei=np.zeros((len(r),3))
    Hamitonian_Zee=np.ones((len(r),3))
    Hamitonian_Zee=-Bi*Hamitonian_Zee
    for ii in np.arange(len(r)):
        for jj in np.arange(NN[ii].astype('int32')):
            Hamitonian_Hei[ii,:]+=(-J1*spin[ii,:]*spin[jj,:])
    
    Hamitonian=Hamitonian_Hei+Hamitonian_Zee                
    return Hamitonian


def find_Beff(spin):
    H=Hamitonian(spin)
    B_eff=np.divide(H,spin,np.zeros_like(H,dtype=np.float64),where=spin!=0)
    damping=find_damping(spin)
    B_eff_damp=np.zeros((len(r),3))
    for ii in np.arange(len(r)):
        for jj in np.arange(len(r)):
            B_eff_damp[ii,:]+=(damping[ii,jj]*derivate(spin[jj,:],B_eff[jj,:]))/np.sqrt(np.sum(np.square(spin[jj,:])))
    B_eff_total=B_eff+B_eff_damp
    return B_eff_total
def find_damping(spin):
    damping=np.identity(len(r))*damping_onsite
    H=Hamitonian(spin)
    for ii in np.arange(len(r)):
            for jj in np.arange(NN[ii].astype('int32')):
                damping[ii,NL1[ii,jj].astype('int32')]=damping_offsite
    return damping

def giromagneticRatio():
    return -1 

def derivate(S,H):
    return giromagneticRatio()*(np.cross(S.T,H.T))


def update_spin(spin, n):
    B_eff=find_Beff(spin)
    result = (derivate(spin[n,:], B_eff[n,:]))
    spin[n,0] = spin[n,0] + dt*result[0]
    spin[n,1] = spin[n,1] + dt*result[1]
    spin[n,2] = spin[n,2] + dt*result[2]
    return spin

def plotMagnetization(mx, my, mz):
	fig, ax = plt.subplots(figsize=(10,10))
	plt.plot(mx, label = 'mx')
	plt.plot(my, label = 'my')
	plt.plot(mz, label = 'mz')
	plt.title('Magnetization')
	plt.xlabel('step',fontdict={'weight': 'normal', 'size': 13})
	plt.ylabel('M',fontdict={'weight': 'normal', 'size': 13})
	plt.legend(loc=2)
  plt.show()

def plotPositions():
    fig = plt.subplots(figsize=(10,10))
    ax = fig.gca(projection='3d')
    x,y,z= np.array((r[:,0],r[:,1], r[:,2]))
    u=spinLattice[:,0]
    v=spinLattice[:,1]
    w=spinLattice[:,2]
    plt.xlim(-2, 5) 
    plt.ylim(-2, 5)
    ax.quiver(x, y,z, u, v, w, color='b')
    plt.show()
        
mx = np.zeros((len(r),step))
my = np.zeros((len(r),step))
mz = np.zeros((len(r),step))
spinLattice=gen_rand_spinlattice()
# plotPositions()
for t in range(step):
 	for s in range(len(r)):
  	 	spin = update_spin(spinLattice, s)
  	 	mx[s][t] = spin[s,0]
  	 	my[s][t] = spin[s,1]
  	 	mz[s][t] = spin[s,2]
spinIndex = 0 
Mx=mx.mean(axis=0)
My=my.mean(axis=0)
Mz=mz.mean(axis=0)
plotMagnetization(mx[spinIndex],my[spinIndex],mz[spinIndex])

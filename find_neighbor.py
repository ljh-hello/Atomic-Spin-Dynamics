# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:25:05 2021

@author: 59778
"""


# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:32:23 2021

@author: 59778
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 22:49:47 2021

@author: luzhiwei
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 22:01:10 2021

@author: 59778
"""



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from expand_cell import expand_cell

def read_poscar(pos_file):
    try:
        latt_vec=np.zeros([3,3])    
        poscar=open(pos_file)
        line=poscar.readline()
        sys_name=line
        latt= float(poscar.readline())
        for i in range(3):
        	line = poscar.readline().split()
        	for j in range(3):
        		latt_vec[i,j] = latt*float(line[j])
        atom_ele=poscar.readline().split()
        atom_num_old=poscar.readline().split()
        atom_num=[]
        for n in atom_num_old:
            atom_num.append(int(n)) 
        atom=atom_ele+atom_num
        atom_sum=(np.array(atom_num)).sum()
        line2=poscar.readline()
        if line2.strip() == 'Selective dynamics':
            line2=poscar.readline()
        S=np.zeros([atom_sum,3])
        for ii in range(atom_sum):
        	S1 = poscar.readline().split()
        	for jj in range(3):
        		S[ii,jj] =np.float(S1[jj])
        if line2.strip()=='D' or 'Direct':
            S=np.dot(S,latt_vec)
        else:
            S=S  
    except IOError:
    	print ('POSCAR is not exist or file open failed! \n')
    poscar.close()
    return latt,latt_vec,atom,S,sys_name,line2  
 
def expand_cell(pos_file,Nx,Ny,Nz):  
    n=0
    latt,latt_vec,atom,S,sys_nam,line2=read_poscar(pos_file)
    latt_s=np.sqrt(np.sum(np.square(latt_vec*np.array((Nx,Ny,Nz)).T),axis=1))
    r=np.zeros((len(S)*Nx*Ny*Nz,3))
    for ii in np.arange(Nx):
        for jj in np.arange(Ny):
            for kk in np.arange(Nz):
                for mm in np.arange(len(S)):
                    r[n,:]=np.dot(np.array([ii,jj,kk]).reshape(1,3),latt_vec)+S[mm,:]
                    n=n+1
    return r,latt_s

def plot_atom(r):
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(r[:,0], r[:,1], r[:,2],c = 'r', marker = 'o',s=100)
    ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red'})
    ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    plt.show()


def find_neighbor_pbc(pos_file,Nx,Ny,Nz):
    r,latt_s=expand_cell(pos_file,Nx,Ny,Nz)
    NN=np.zeros((len(r),1))
    NL1=(np.zeros((len(r),len(r))))
    N_rc=np.zeros((len(r)-1,1))
    for ii in np.arange(1,len(r)):
        dis1=np.sqrt(np.sum(np.square((r[0,:]-r[ii,:]))))
        N_rc[ii-1,:]=dis1.astype('float')
    N_rc2=np.sort(N_rc,axis=0)
    rc=N_rc2[0]+0.1
    # print(N_rc2)
    for n1 in np.arange(0,len(r)):
        for n2 in np.arange(0,len(r)):
            r12=np.sqrt(np.sum(np.square((r[n1,:]-r[n2,:])-np.round((r[n1,:]-r[n2,:])/latt_s)*latt_s)))
            if r12<rc and r12!=0:
                NN[n1]=NN[n1]+1
                NL1[n1,NN[n1].astype('int64')-1]=n2
    # print(np.sqrt(np.sum(np.square((r[0,:]-r[11,:])-np.round((r[0,:]-r[11,:])/latt_s)*latt_s))))
    return NN,NL1



def find_neighbor_npbc(pos_file,Nx,Ny,Nz):
    r,latt_s=expand_cell(pos_file,Nx,Ny,Nz)
    NN=np.zeros((len(r),1))
    NL1=(np.zeros((len(r),len(r))))
    N_rc=np.zeros((len(r)-1,1))
    for ii in np.arange(1,len(r)):
        dis1=np.sqrt(np.sum(np.square((r[0,:]-r[ii,:]))))
        N_rc[ii-1,:]=dis1.astype('float')
    N_rc2=np.sort(N_rc,axis=0)
    rc=N_rc2[0]+0.1
    for n1 in np.arange(0,len(r)):
        for n2 in np.arange(0,len(r)):
            r12=np.sqrt(np.sum(np.square(r[n1,:]-r[n2,:])))
            # print(dis)
            if r12<rc and r12!=0:
                NN[n1]=NN[n1]+1
                NL1[n1,NN[n1].astype('int64')-1]=n2
    
    return NN,NL1



























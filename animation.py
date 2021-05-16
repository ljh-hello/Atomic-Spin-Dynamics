import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from pathlib import Path

# data = np.loadtxt("/Users/luzhiwei/Research/fortran_study/task_nonlocal_damping/fortran_new_danny/fortran/test/Dimer/damping_type_different/damp_type_1/moment.dat")

data = np.loadtxt("/Users/luzhiwei/Research/fortran_study/task_nonlocal_damping/fortran_new_danny/fortran/test/Dimer/damping_type_different/new_No_field/non_1.0/damp_type_5/moment.dat")
fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))
X=np.array([0,0.3])
Y=np.array([0,0])
Z=np.array([0,0])
U = data[0:2,3]
V = data[0:2,4]
W = data[0:2,5]
plt.xlim(-0.2,0.3) 
plt.ylim(-1,1)
Q=ax.quiver(X,Y,Z,U,V,W,length=0.05,color='b', normalize=True,arrow_length_ratio=0.3,linewidths=2)


def update_quiver(num):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    global Q
    Q.remove()
    
    U = data[num:num+2,3]
    V = data[num:num+2,4]
    W =  data[num:num+2,5]
    Q=ax.quiver(X,Y,Z,U,V,W,length=0.05,color='b', normalize=True,arrow_length_ratio=0.3,linewidths=2)
    return Q,


anim = animation.FuncAnimation(fig, update_quiver, frames=[i*2 for i in range(0,int(0.5*data.shape[0]))],interval=1, blit=False)
fig.tight_layout()

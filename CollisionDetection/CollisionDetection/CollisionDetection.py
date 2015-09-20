# 2015.09.20, LI Yunsheng

import numpy as np
from scipy.fftpack import fft2, ifft2
import matplotlib.pyplot as plt


class Grid():
    # M - column, N - row must be even
    def __init__(self,dM,M,dN,N,mat):
        self.M = M
        self.N = N
        self.dM = dM
        self.dN = dN
        self.mat = mat
    @property
    def width(self):
        return M*dM
    @property
    def height(self):
        return N*dN


class Circle():
    def __init__(self, r):
        #self.x = x
        #self.y = y
        self.r = r
    def mesh(self, grid_map):
        # return NXM array
        R = np.zeros((grid_map.N, grid_map.M)) 
        j = np.floor(-self.r/grid_map.dN + 0.5)
        y = (j+0.5)*grid_map.dN
        while j < 0:
            x = np.sqrt(self.r**2 - y**2)
            i = np.floor(-x/grid_map.dM + 0.5)
            i_list = np.linspace(i, -i, -2*i+1)
            i_list = np.where(i_list<0, i_list+grid_map.M, i_list)
            jj=j+grid_map.N
            for ii in i_list:
                R[ii,jj] = 1
            j += 1
            y += grid_map.dN
        j = 0
        y = -grid_map.dN/2
        while y < self.r:
            x = np.sqrt(self.r**2 - y**2)
            i = np.floor(-x/grid_map.dM + 0.5)
            i_list = np.linspace(i, -i, -2*i+1)
            i_list = np.where(i_list<0, i_list+grid_map.N, i_list)
            for ii in i_list:
                R[ii,j] = 1
            j += 1
            y += grid_map.dN
        return R

class Veh_Cfg():
    def __init__(self, x,y,t,l1,l2,w):
        self.x = x
        self.y = y
        self.t = t
        self.l1 = l1
        self.l2 = l2
        self.w = w


'''
disk = Circle(10)
grdmap = Grid(0.1,500,0.1,500,np.zeros((500,500)))
data = disk.mesh(grdmap)
plt.imshow(data, cmap=plt.cm.gray_r)
plt.show()
'''

N=400 # map size
delta =0.25 # incremental size
eps = 0.1 # err

o = np.zeros((N,N)) # obstacles
o[0, :] = 1
o[N-1, :] = 1
o[:, 0] = 1
o[:, N-1] = 1
o[170:230, 170:230] = 1
obstacles = Grid(delta,N,delta,N,o)

disk = Circle(10)
robot = disk.mesh(obstacles)

Obst = fft2(obstacles.mat)
Robot = fft2(robot)
Collision = Obst * Robot
collision = np.real(ifft2(Collision))
collision = np.where(collision > eps, 0.5, 0) + obstacles.mat
collision = np.where(collision > 1, 1, collision)
plt.imshow(collision, cmap=plt.cm.gray_r)
plt.show()
# r = np.zeros((N, N)) # robot
# r[0:5, 0:5] = 1
# r[N-4:N, 0:5]= 1
# r[0:5, N-4:N] = 1
# r[N-4:N, N-4:N] = 1
# O = fft2(o)
# R = fft2(r)
# C = O*R # convolution
# c = np.real(ifft2(C))
# c = np.where(c > eps, 1, 0)
# plt.imshow(c, cmap=plt.cm.gray_r)
# plt.colorbar()
# plt.show()

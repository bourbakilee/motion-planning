# 2015.09.20, LI Yunsheng

import numpy as np
from scipy.fftpack import fft2, ifft2
import matplotlib.pyplot as plt


class Grid():
    # M - column, N - row must be even
    def __init__(self,dM,dN,data):
        self.M = data.shape[1]
        self.N = data.shape[0]
        self.dM = dM
        self.dN = dN
        self.data = data

    @property
    def width(self):
        return M*dM

    @property
    def height(self):
        return N*dN


class Circle():
    def __init__(self, r):
        #self.x = 0
        #self.y = 0
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
            #i_list = np.where(i_list<0, i_list+grid_map.M, i_list)
            #jj=j+grid_map.N
            for ii in i_list:
                R[ii,j] = 1
            j += 1
            y += grid_map.dN
        j = 0
        y = -grid_map.dN/2
        while y < self.r:
            x = np.sqrt(self.r**2 - y**2)
            i = np.floor(-x/grid_map.dM + 0.5)
            i_list = np.linspace(i, -i, -2*i+1)
            #i_list = np.where(i_list<0, i_list+grid_map.N, i_list)
            for ii in i_list:
                R[ii,j] = 1
            j += 1
            y += grid_map.dN
        return R

    def moveto(self,disk_mesh,x,y):
        # disk_mesh - NXM array
        # x,y > 0
        R = np.zeros(disk_mesh.shape)
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                R[i,j] = disk_mesh[i-x, j-y] 
        return R


class Veh_Cfg():
    # (x,y)-center of rear axis
    def __init__(self, x,y,t,l1,l2,w):
        self.x = x
        self.y = y
        self.t = t
        self.l1 = l1
        self.l2 = l2
        self.length = l1+l2
        self.width = w
        self.r = np.sqrt(self.length**2/9 + self.width**2/4) # radius of circles, which cover the vehicle
        self.d = 2*self.length/3

    def centers_of_circles(self,grid_map):
        c = np.zeros((3,2))
        direction = np.array([np.sin(self.t),np.cos(self.t)])
        c[1] = np.array([self.x,self.y]) + (self.l1 - self.l2)/2 * direction
        c[0] = c[1] - self.d * direction
        c[2] = c[1] + self.d * direction
        return np.floor(c/np.array([[grid_map.dM, grid_map.dN]]))

    def cost(self, centers, costmap):
        return np.max([costmap[tuple(c)] for c in centers])



N=1000 # map size
delta =0.1 # incremental distance
eps = 0.1 # numerical err

obstacles = np.zeros((N,N)) # obstacles
obstacles[0, :] = 1
obstacles[N-1, :] = 1
obstacles[:, 0] = 1
obstacles[:, N-1] = 1
obstacles[400:600, 400:600] = 1
workspace = Grid(delta,delta,obstacles)

veh = Veh_Cfg(25,25,np.pi/4,4,1,2)
disk = Circle(veh.r)
disk_mesh = disk.mesh(workspace)
centers = veh.centers_of_circles(workspace)
veh_mesh = np.zeros((N,N))
veh_mesh += disk.moveto(disk_mesh,centers[0,0],centers[0,1])
veh_mesh += disk.moveto(disk_mesh,centers[1,0],centers[1,1])
veh_mesh += disk.moveto(disk_mesh,centers[2,0],centers[2,1])
veh_mesh = np.where(veh_mesh>0, 0.8, 0)
# 还是Emacs好用，呵呵
Obstacles = fft2(obstacles)
Robot = fft2(disk_mesh)
CostMap = Obstacles * Robot
costmap = np.real(ifft2(CostMap))
costmap = np.where(costmap > eps, 0.5, 0) + obstacles + veh_mesh
costmap = np.where(costmap > 1, 1, costmap)

cost = veh.cost(centers,costmap)

print(cost)
x=np.linspace(1,99,990)
y=20+10*np.sin(x)
plt.imshow(costmap, cmap=plt.cm.gray_r, origin="lower", extent=(0,100,0,100))
plt.plot(x,y)
# plt.rc('figure', figsize=(1000,1000))
plt.show()

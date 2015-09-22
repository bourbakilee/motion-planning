# 2015.09.21, LI Yunsheng

import numpy as np
import spiral
from scipy.fftpack import fft2, ifft2
import matplotlib.pyplot as plt

class Road():
    def __init__(self, length, lane_num=2, lane_width=3.75, center_line=None, center_line_fun=None, q=None, ref_grid_width=0.6, ref_grid_length=1, ref_delta_s=0.1):
        # center_line: Nx4 array, [[x,y,theta,k],...]
        self.length = length
        self.lane_num = lane_num
        self.lane_width = lane_width
        self.width = lane_num*lane_width
        self.grid_num_per_lane = int(2 * np.ceil(np.ceil(lane_width / ref_grid_width) / 2)) # lateral direction
        self.grid_num_lateral = int(self.grid_num_per_lane*lane_num)
        self.grid_num_longitudinal = int(np.ceil(length/ref_grid_length))
        self.grid_width = lane_width / self.grid_num_per_lane
        self.grid_length = length/self.grid_num_longitudinal
        self.N = int(np.ceil(length/ref_delta_s)) # number of discretized points of center line
        self.delta_s = length/(self.N-1)
        if center_line_fun is not None:
            self.center_line = spiral.spiral_calc(center_line_fun, length, q=q, ref_delta_s=ref_delta_s) # q0=(0,0,0)
        elif center_line is not None:
            self.center_line = center_line # do not use this
        else:
            line = np.zeros((self.N,4))
            line[:,0] = np.linspace(0,length,self.N)
            if q is not None:
                sin_x = np.sin(q[2])*line[:,0]
                cos_x = np.cos(q[2])*line[:,0]
                sin_y = np.sin(q[2])*line[:,1]
                cos_y = np.cos(q[2])*line[:,1]
                line[:,0] = q[0] + cos_x - sin_y
                line[:,1] = q[1] + sin_x + cos_y
                line[:,2] = np.mod((line[:,2] + q[2]), 2*np.pi)
            self.center_line = line
        boundaries = np.zeros((self.N, 2*(lane_num+1)))
        lateral_bias = np.linspace(-lane_num*lane_width/2, lane_num*lane_width/2, lane_num+1)
        for i in range(lane_num+1):
            boundaries[:,(2*i):(2*i+2)] = self.lateral_biasing_line(lateral_bias[i])
        self.boundary_lines = boundaries

    def lateral_biasing_line(self,lateral_bias):
        # lateral_bias \in (-self.width/2, self.width/2)
        # return Nx2 array (x,y)
        return self.center_line[:,0:2] + lateral_bias*np.array([-np.sin(self.center_line[:,2]), np.cos(self.center_line[:,2])]).T

    def longitudinal_biasing_line(self,longitudinal_bias):
        # longitudinal_bias \in (0,self.length)
        # return Mx2 array (x,y)
        lateral_bias = np.linspace(-self.width/2, self.width/2,np.ceil(self.width/0.1))
        base_station = self.sl2xy(longitudinal_bias, 0) #(x,y,theta,k)
        return base_station[0:2] + np.array([-np.sin(base_station[2])*lateral_bias, np.cos(base_station[2])*lateral_bias]).T

    def sl2xy(self,s,l):
        # s - longitudinal bias, l - lateral bias
        if np.abs(s)<1.e-4:
            station = self.center_line[0]
        elif np.abs(s-self.length)<1.e-4:
            station = self.center_line[-1]
        else:
            k = np.floor(s/self.length * self.N)
            w1 = s/self.delta_s - k
            w2 = 1 - w1
            station = self.center_line[k]*w2 + self.center_line[k+1]*w1
        return np.array([station[0]-l*np.sin(station[2]), station[1]+l*np.cos(station[2]), station[2], (station[3]**-1-l)**-1])

    def ij2xy(self,i,j):
        # j - lateral offset \in [-self.grid_num_lateral/2, self.grid_num_lateral.2]
        # i - longitudinal offset \in [0,self.grid_num_longitudinal]
        return self.sl2xy(i*self.grid_length, j*self.grid_width)

class Vehicle():
    def __init__(self, cfg=np.array([0.,0.,0.]),length=5., width=2.,wheelbase=3.):
        #trajectory:Nx4 array - (t,x,y,theta)
        self.cfg = cfg # Confiduration
        self.center_of_rear_axle = cfg[0:2]
        self.head = cfg[2]
        self.length = length
        self.width = width
        self.wheelbase = wheelbase
        c, s = np.cos(cfg[2]), np.sin(cfg[2])
        self.geometric_center = cfg[0:2]+wheelbase/2*np.array([c,s])
        self.vertex = self.geometric_center + 0.5*np.array([
            [c*length-s*width, s*length+c*width],
            [-c*length-s*width,-s*length+c*width],
            [-c*length+s*width,-s*length-c*width],
            [c*length+s*width, s*length-c*width]
            ])
        # self.trajectory = trajectory if trajectory is not None else np.array([-1,0,0,0])

    @property
    def covering_disk_radius(self):
        return  np.sqrt(self.length**2/9. + self.width**2/4.)

    @property
    def covering_disk_centers(self):
        distance = 2.*self.length/3.
        direction = np.array([np.sin(self.head),np.cos(self.head)])
        centers = np.zeros((3,2))
        centers[1] = self.geometric_center
        centers[0] = centers[1] - distance * direction
        centers[2] = centers[1] + distance * direction
        return centers

    def moveto(self,q):
        return Vehicle(cfg=q,length=self.length, width=self.width,wheelbase=self.wheelbase)


class Workspace():
    def __init__(self, resolution=0.25, M=320, N=400, static_obst=None, moving_obst=None, road=None,vehicle=Vehicle(cfg=np.array([3.,10.,0.]))):
        # M - row number of workspace
        # N - column number of workspace
        # static_obst - MxN arrary of {0 - free,1 - occupied}
        # moving_obst - list of vehicle shape, optional
        # road - optional
        # vehicle - planning object
        self.resolution = resolution
        self.M = M
        self.N = N
        self.static_obst = np.zeros((M,N)) if static_obst is None else static_obst
        self.moving_obst = moving_obst
        self.road = road
        self.vehicle = vehicle
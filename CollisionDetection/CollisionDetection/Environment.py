# 2015.11.16, LI Yunsheng

from math import *
import cv2
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import fsolve

class Vehicle():
    def __init__(self, length=4., width=1.6, wheelbase=2.4, trajectory=np.array([[-1.,-1.,50.125,50.125,0.,0.,0.,0.,0.,0.]])):
        # trajectory: N X 10 array - [ [t, s, x, y, theta, k, dk, v, a, j] ]
        # state: 1 X 6 vector - [t, x, y, theta, k, v]
        self.trajectory = trajectory
        self.state = np.array([trajectory[0,0], trajectory[0,2], trajectory[0,3], trajectory[0,4], trajectory[0,5], trajectory[0,7]])
        self.traj_fun = None if trajectory[0,0] < 0. else interp1d(self.trajectory[:,0],self.trajectory[:,1:].T,kind='linear')
        # self.traj_fun2 = None if trajectory[0,1] < 0. else interp1d(self.trajectory[:,1],np.array([self.trajectory[:,0],self.trajectory[:,2],self.trajectory[:,3],self.trajectory[:,4],self.trajectory[:,5],self.trajectory[:,6],self.trajectory[:,7],self.trajectory[:,8],self.trajectory[:,9]]).T)
        # geometric parameters
        self.length = length
        self.width = width
        self.wheelbase = wheelbase
        #
        self.center_of_rear_axle = self.state[1:3]
        self.heading = self.state[3]
        c, s = np.cos(self.heading), np.sin(self.heading)
        self.geometric_center = self.center_of_rear_axle + wheelbase/2*np.array([c,s])
        self.vertex = self.geometric_center + 0.5*np.array([
            [-c*length-s*width, -s*length+c*width],
            [-c*length+s*width,-s*length-c*width],
            [c*(length-0.3*width)+s*width, s*(length-0.3*width)-c*width],
            [c*length+s*0.7*width,s*length-c*0.7*width],
            [c*length-s*0.7*width, s*length+c*0.7*width],
            [c*(length-0.3*width)-s*width, s*(length-0.3*width)+c*width]
            ])


    def update(self,time):
        if self.traj_fun is not None and time is not None and time <= self.trajectory[-1,0]:
            tmp = self.traj_fun(time) # [s, x, y, theta, k, dk, v, a, j]
            self.state = np.array([time, tmp[1], tmp[2], tmp[3], tmp[4], tmp[6]])
            self.center_of_rear_axle = self.state[1:3]
            self.heading = self.state[3]
            c, s = np.cos(self.heading), np.sin(self.heading)
            self.geometric_center = self.center_of_rear_axle + wheelbase/2*np.array([c,s])
            self.vertex = self.geometric_center + 0.5*np.array([
                [-c*length-s*width, -s*length+c*width],
                [-c*length+s*width,-s*length-c*width],
                [c*(length-0.3*width)+s*width, s*(length-0.3*width)-c*width],
                [c*length+s*0.7*width,s*length-c*0.7*width],
                [c*length-s*0.7*width, s*length+c*0.7*width],
                [c*(length-0.3*width)-s*width, s*(length-0.3*width)+c*width]
                ])
        # else:
        #     print("invalid time:", time)


    def recover(self):
        trajectory = self.trajectory
        length = self.length
        width = self.width
        wheelbase = self.wheelbase
        self.state = np.array([trajectory[0,0], trajectory[0,2], trajectory[0,3], trajectory[0,4], trajectory[0,5], trajectory[0,7]])
        self.center_of_rear_axle = self.state[1:3]
        self.heading = self.state[3]
        c, s = np.cos(self.heading), np.sin(self.heading)
        self.geometric_center = self.center_of_rear_axle + wheelbase/2*np.array([c,s])
        self.vertex = self.geometric_center + 0.5*np.array([
            [-c*length-s*width, -s*length+c*width],
            [-c*length+s*width,-s*length-c*width],
            [c*(length-0.3*width)+s*width, s*(length-0.3*width)-c*width],
            [c*length+s*0.7*width,s*length-c*0.7*width],
            [c*length-s*0.7*width, s*length+c*0.7*width],
            [c*(length-0.3*width)-s*width, s*(length-0.3*width)+c*width]
            ])

    #
    def covering_disk_radius(self):
        return np.sqrt(self.length**2/9. + self.width**2/4.)

    #
    def covering_disk_centers(self):
        distance = 2.*self.length/3.
        direction = np.array([np.sin(self.heading),np.cos(self.heading)])
        centers = np.zeros((3,2))
        centers[1] = self.geometric_center
        centers[0] = centers[1] - distance * direction
        centers[2] = centers[1] + distance * direction
        return centers



class Road():
    def __init__(self, center_line, lane_width=3.75, ref_grid_width=0.6, ref_grid_length=1.2):
        # center_line: Nx5 array, [[s, x, y, theta, k],...]
        # lane_number = 3
        self.length = center_line[-1,0]
        self.lane_width = lane_width
        self.width = lane_width * 3.
        self.center_line = center_line 
        self.center_line_fun = interp1d(center_line[:,0], center_line[:,1:].T, kind='linear') # return [x,y,theta,k]
        self.grid_num_per_lane = 2 * ceil(ceil(lane_width / ref_grid_width) / 2) # lateral direction
        self.grid_num_lateral = self.grid_num_per_lane*3 # 横向网格数目
        self.grid_num_longitudinal = ceil(self.length/ref_grid_length) # 纵向网格数目
        self.grid_width = lane_width / self.grid_num_per_lane
        self.grid_length = self.length/self.grid_num_longitudinal
        #
        self.lateral_biases = np.linspace(-self.width/2., self.width/2., self.grid_num_lateral+1)
        self.longitudinal_biases = np.linspace(0., self.length, self.grid_num_longitudinal+1)
        # 计算纵横向网格线，方便绘图
        longitudinal_lines = np.zeros((self.center_line.shape[0], 2*(self.grid_num_lateral+1))) #[[x,y],...]
        for i in range(self.grid_num_lateral+1):
            longitudinal_lines[:,(2*i):(2*i+2)] = self.lateral_biasing_line(self.lateral_biases[i])
        self.longitudinal_lines = longitudinal_lines

        lateral_lines = np.zeros((ceil(self.width/0.1), 2*(self.grid_num_longitudinal+1))) #[[x,y],...]
        for i in range(self.grid_num_longitudinal+1):
            lateral_lines[:,(2*i):(2*i+2)] = self.longitudinal_biasing_line(self.longitudinal_biases[i])
        self.lateral_lines = lateral_lines
        #


    def lateral_biasing_line(self,lateral_bias):
        # lateral_bias \in (-self.width/2, self.width/2)
        # return Nx2 array (x,y): x=x0-l*sin(theta), y=y0+l*cos(theta)
        return self.center_line[:,1:3] + lateral_bias*np.array([-np.sin(self.center_line[:,3]), np.cos(self.center_line[:,3])]).T


    def longitudinal_biasing_line(self,longitudinal_bias):
        # longitudinal_bias \in (0,self.length)
        # return Mx2 array (x,y)
        lateral_bias = np.linspace(-self.width/2., self.width/2.,np.ceil(self.width/0.1))
        base_station = self.sl2xy(longitudinal_bias, 0) #(x,y,theta,k)
        return base_station[0:2] + np.array([-np.sin(base_station[2])*lateral_bias, np.cos(base_station[2])*lateral_bias]).T


    def sl2xy(self, s, l):
        # s - longitudinal bias, l - lateral bias
        # return [x,y,theta,k]
        if 0<=s<=self.length:
            station = self.center_line_fun(s)
            return np.array([station[0]-l*np.sin(station[2]), station[1]+l*np.cos(station[2]), station[2], (station[3]**-1-l)**-1])
        else:
            return None


    def ij2xy(self,i,j):
        # return [x,y,theta,k]
        if abs(j) <= self.grid_num_lateral//2 and 0<=i<=self.grid_num_longitudinal:
            return self.sl2xy(self.longitudinal_biases[i], self.lateral_biases[j+self.grid_num_lateral//2])
        else:
            return None


    def __xys(self,x,y,s):
        tmp = self.center_line_fun(s)
        return (x-tmp[0])*np.cos(tmp[2]) + (y-tmp[1])*np.sin(tmp[2])
    

    def xy2sl(self, x, y):
        f = lambda s: self.__xys(x,y,s)
        s0 = fsolve(f, self.length/2)
        tmp = self.center_line_fun(s0)
        if abs(tmp[2]) < 1.e-4:
            l0 = (y-tmp[1])/np.cos(tmp[2])
        else:
            l0 = (tmp[0]-x)/np.sin(tmp[2])
        return np.array([s0[0], l0[0]])



class Workspace():
    def __init__(self,base=np.zeros((401,401)), resolution=0.25, vehicle=Vehicle(), road=None, lane_costs=None, static_obsts=None, moving_obsts=None):
        self.base = base
        self.row = base.shape[0]
        self.column = base.shape[1]
        self.resolution = resolution
        self.static_obsts = static_obsts # list of static vehicles
        self.moving_obsts = moving_obsts # list of moving vehicles
        self.road = road
        self.lane_grids = self.grids_of_lanes(self.road)
        self.lane_costs = lane_costs
        self.lane_map = self.__lane_map()
        self.vehicle = vehicle
        #
        self.disk = self.disk_filter()
        #
        self.time = 0. if self.moving_obsts else None
        self.static_map = self.__static_map()
        self.env_map = self.__env_map()
        self.collision_map = self.__collision_map()
        self.cost_map = self.__cost_map()


    def update(self, time):
        if time is not None and time >= 0:
            self.time = time
            for i in range(len(self.moving_obsts)):
                self.moving_obsts[i].update(time)
            self.env_map = self.__env_map()
            self.collision_map = self.__collision_map()
            self.cost_map = self.__cost_map()


    def set_lane_costs(self, lane_costs):
        # lane_costs: [v_right, v_center, v_left]
        self.lane_costs = lane_costs
        self.lane_map = self.__lane_map()
        self.cost_map = self.__cost_map()


    def disk_filter(self, r=None):
        # r <= 10m, O(10.125, 10.125)
        if r is None:
            r = self.vehicle.covering_disk_radius()
        R = np.zeros((81,81))
        k0 = ceil(r/self.resolution - 0.5)
        # N = 2 * k0 + 1 # 滤波器是NxN方阵，N为奇数
        # 中间一行
        for j in range(40-k0, 41+k0):
            R[40,j] = 1.
        #向上\下
        for i in range(0, k0):
            rr = sqrt( r**2 - ((i+0.5)*self.resolution)**2 )
            k = ceil(rr/self.resolution - 0.5)
            for j in range(40-k, 41+k):
                R[41+i, j] = 1.
                R[39-i, j] = 1.
        F = R[(40-k0):(41+k0), (40-k0):(41+k0)]
        return F


    def vehicle_filter(self,theta=0.):
        veh = Vehicle(trajectory=np.array([[-1.,-1.,50.125,50.125,theta,0.,0.,0.,0.,0.]]))
        #R = np.zeros((81,81))
        r = sqrt(veh.width**2 + (veh.length + veh.wheelbase)**2) / 2
        k0 = ceil(r/self.resolution - 0.5)
        #
        # grids_list = []
        # for i in range(len(veh.vertex)):
        #     j = (i+1) % (len(veh.vertex))
        #     #
        #     j_m1 = floor(veh.vertex[i,0]/self.resolution)
        #     j_m2 = floor(veh.vertex[j,0]/self.resolution)
        #     if j_m1 == j_m2:
        #         grids = [(j_m1,floor(veh.vertex[i,1]/self.resolution)), (j_m2,floor(veh.vertex[j,1]/self.resolution))]
        #     else:
        #         f = interp1d([veh.vertex[i,0], veh.vertex[j,0]], [veh.vertex[i,1], veh.vertex[j,1]])
        #         if veh.vertex[i,0] < veh.vertex[j,0]:
        #             grids = self.grids_occupied_by_line(f,veh.vertex[i,0], veh.vertex[j,0])
        #         else:
        #             grids = self.grids_occupied_by_line(f,veh.vertex[j,0], veh.vertex[i,0])
        #     grids_list.append(grids)
        F = self.grids_occupied_by_polygon(veh.vertex)
        return F[(200-k0):(201+k0),(200-k0):(201+k0)]


    def grids_occupied_by_line(self, f, x1, x2):
        # 0 <= x1 < x2 <= 100
        j_min, j_max = floor(x1/self.resolution), floor(x2/self.resolution) #列
        # f 单调 : 0 <= f(x1) , f(x2) <= 100
        i_min, i_max = floor(f(x1)/self.resolution), floor(f(x2)/self.resolution) #行
        #
        if x1 > x2:
            j_min, j_max = j_max, j_min
            i_min, i_max = i_max, i_min
        #
        grids = [(j_min, i_min), (j_max, i_max)]
        for j in range(j_min+1, j_max+1):
            i = floor(f(j*self.resolution)/self.resolution)
            grids.append((j-1,i))
            grids.append((j,i))
        return grids #[(列，行)]


    def grids_encircled_by_lines(self,grids_list):
        R = np.zeros((401,401))
        grids = []
        for i in range(len(grids_list)):
            grids.extend(grids_list[i])
        grids.sort()
        i = 0
        j = 1
        while i < len(grids):
            while grids[j][0] == grids[i][0]:
                j += 1
                if j == len(grids):
                    break
            for k in range(grids[i][1], grids[j-1][1]+1):
                R[k,grids[i][0]] = 1.
            if j == len(grids):
                break
            else:
                i = j
                j += 1
        return R


    def grids_occupied_by_polygon(self,plg):
        # plg : List of vertex [(x1,y1),(x2,y2),...,(xn,yn)]
        grids_list=[]
        for i in range(len(plg)):
            j = (i+1) % (len(plg))
            j_m1 = floor(plg[i][0] / self.resolution)
            j_m2 = floor(plg[j][0] / self.resolution)
            if j_m1 == j_m2:
                grids = [(j_m1,floor(plg[i][1]/self.resolution)), (j_m2,floor(plg[j][1]/self.resolution))]
            else:
                f = interp1d([plg[i][0], plg[j][0]], [plg[i][1], plg[j][1]])
                # if plg[i][0] < plg[j][0]:
                grids = self.grids_occupied_by_line(f,plg[i][0], plg[j][0])
                # else:
                #     grids = self.grids_occupied_by_line(f,plg[j][0], plg[i][0])
            grids_list.append(grids)
        return self.grids_encircled_by_lines(grids_list)


    def grids_of_lanes(self, road):
        # return: list of matrix map
        # if road is None:
        #     return None
        # else:
        ll = road.longitudinal_lines
        nl = road.grid_num_per_lane
        # print(nl)
        map_list = []
        f1 = interp1d(ll[:,0], ll[:,1])
        for i in range(3):
            grids_list = []
            f2 = interp1d(ll[:,2*(i+1)*nl], ll[:,1+2*(i+1)*nl])
            g = interp1d([ll[0,2*i*nl], ll[0,2*(i+1)*nl]], [ll[0, 1+2*i*nl], ll[0, 1+2*(i+1)*nl]])
            h = interp1d([ll[-1,2*i*nl], ll[-1,2*(i+1)*nl]], [ll[-1, 1+2*i*nl], ll[-1, 1+2*(i+1)*nl]])
            # print(2*i*nl, 2*(i+1)*nl, 1+2*i*nl, 1+2*(i+1)*nl)
            # print(ll[0,2*i*nl], ll[0,2*(i+1)*nl], ll[0, 1+2*i*nl], ll[0, 1+2*(i+1)*nl])
            # print(g)
            # 可以修改grids_occupied_by_line函数，自动检测参数的大小，避免下面这样的处理。
            # if ll[0,2*i*nl] < ll[0,2*(i+1)*nl]:
            grids1 = self.grids_occupied_by_line(g, ll[0,2*i*nl], ll[0,2*(i+1)*nl])
            # else:
            #     grids1 = self.grids_occupied_by_line(g, ll[0,2*(i+1)*nl], ll[0,2*i*nl])
            # if ll[-1,2*i*nl] < ll[-1,2*(i+1)*nl]:
            grids2 = self.grids_occupied_by_line(h, ll[-1,2*i*nl], ll[-1,2*(i+1)*nl])
            # else:
            #     grids2 = self.grids_occupied_by_line(h, ll[-1,2*(i+1)*nl], ll[-1,2*i*nl])
            # print(grids1)
            # print(grids2)
            # 下面默认f1、f2的参数是从小到大的。
            grids_list.append(self.grids_occupied_by_line(f1, ll[0,2*i*nl], ll[-1,2*i*nl]))
            grids_list.append(self.grids_occupied_by_line(f2, ll[0,2*(i+1)*nl], ll[-1,2*(i+1)*nl]))
            grids_list.append(grids1)
            grids_list.append(grids2)
            map_list.append(self.grids_encircled_by_lines(grids_list))
            f1 = f2
        # print(map_list)
        return map_list


    def __static_map(self):
        if not self.static_obsts:
            return self.base
        else:
            R = self.base
            for i in range(len(self.static_obsts)):
                R += self.grids_occupied_by_polygon(self.static_obsts[i].vertex)
            R = np.where(R > 0., 1., 0.)
            return R


    def __env_map(self):
        R = self.static_map
        if self.moving_obsts is not None:
            for i in range(len(self.moving_obsts)):
                R += self.grids_occupied_by_polygon(self.moving_obsts[i].vertex)
            R = np.where(R > 0., 1., 0.)
        return R


    def __collision_map(self, flt=None):
        if flt is None:
            flt = self.disk
        eps = 1.e-6
        Env_map = self.env_map
        Dest = cv2.filter2D(Env_map, -1, flt)
        Dest = np.where(Dest > eps, 1., 0.)
        return Dest


    def __lane_map(self):
        v1 = min(self.lane_costs)
        p1 = self.lane_costs.index(v1)
        v3 = max(self.lane_costs)
        p3 = self.lane_costs.index(v3)
        p2 = [i for i in [0,1,2] if i!=p1 and i!=p3][0]
        v2 = self.lane_costs[p2]
        lm = v1*self.lane_grids[p1] + v2*self.lane_grids[p2]
        lm = np.where(lm>v2, (v1+v2)/2, lm)
        lm += v3*self.lane_grids[p3]
        lm = np.where(lm>v3, (v2+v3)/2, lm)
        return lm


    def __cost_map(self, flt=None, scale=255.):
        if flt is None:
            flt = self.disk
        eps = 1.e-6
        cost_map = self.collision_map
        cost_map += cv2.filter2D(cost_map, -1, flt)
        cost_map = np.where(cost_map<eps, 0., cost_map)
        cost_map = scale * np.where(cost_map>1, 1., cost_map)
        for i in range(3):
            cost_map += self.lane_costs[i] * self.lane_grids[i]
        cost_map = np.where(cost_map>255., 255., cost_map)
        return cost_map

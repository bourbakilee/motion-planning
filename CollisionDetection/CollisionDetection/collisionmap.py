import numpy as np 
from scipy.interpolate import interp1d
from math import *
import matplotlib.pyplot as plt 
from matplotlib.path import Path
import matplotlib.patches as patches
import cv2
from costmap import Vehicle


def grids_occupied_by_line(f,x1,x2):
    # y=f(x)
    #R = np.zeros((400,400))
    xmin, xmax = floor(x1*4), floor(x2*4)
    ymin, ymax = floor(f(x1)*4), floor(f(x2)*4)
    grid_list = [(xmin, ymin),(xmax,ymax)]
    for xi in range(xmin+1,xmax+1):
        yi = floor(f(xi/4)*4)
        grid_list.append((xi-1,yi))
        grid_list.append((xi,yi))
        #R[yi,xi+1]=1
        #R[yi,xi] = 1
    return grid_list


def grids_encircled_by_line_grids(grid_lists):
    grids = []
    for i in range(len(grid_lists)):
        grids.extend(grid_lists[i])
    grids.sort()
    # print(grids)
    R = np.zeros((81,81))
    i = 0
    j = 1
    # print(len(grids))
    while i < len(grids):
        if i == len(grids)-1:
            break
        #print('i=',i)
        while j < len(grids):
            #print('j=',j)
            if (grids[i][0] != grids[j][0] or j == len(grids)-1):
                #print(grids[i],grids[j-1])
                for k in range(grids[i][1], 1+grids[j-1][1]):
                    R[k, grids[i][0]] = 1
                i = j
                j +=1
                break
            else:
                j += 1
    return R 


def grids_occupied_by_vehicle(veh):
    # R = np.zeros((400,400))
    grid_lists=[]
    for i in range(len(veh.vertex)):
        j = (i+1) % (len(veh.vertex))
        if floor(veh.vertex[i,0]*4) == floor(veh.vertex[j,0]*4):
            grid_list = [(floor(veh.vertex[i,0]*4),floor(veh.vertex[i,1]*4)),(floor(veh.vertex[j,0]*4),floor(veh.vertex[j,1]*4)+1)]
        else:
            f = interp1d([veh.vertex[i,0], veh.vertex[j,0]], [veh.vertex[i,1], veh.vertex[j,1]])
            if veh.vertex[i,0] < veh.vertex[j,0]:
                grid_list = grids_occupied_by_line(f, veh.vertex[i,0], veh.vertex[j,0])
            else:
                grid_list = grids_occupied_by_line(f, veh.vertex[j,0], veh.vertex[i,0])
        grid_lists.append(grid_list)
    return grids_encircled_by_line_grids(grid_lists)


def grids_occupied_by_polygon(plg):
    # plg : List of vertex [(x1,y1),(x2,y2),...,(xn,yn)]
    grid_lists=[]
    for i in range(len(plg)):
        j = (i+1) % (len(plg))
        if np.abs(plg[i][0] - plg[j][0]) < 1.e-8:
            grid_list = [(floor(plg[i][0]*4),floor(plg[i][1]*4)),(floor(plg[j][0]*4),floor(plg[j][1]*4))]
        else:
            f = interp1d([plg[i][0], plg[j][0]], [plg[i][1], plg[j][1]])
            if plg[i][0] < plg[j][0]:
                grid_list = grids_occupied_by_line(f, plg[i][0], plg[j][0])
            else:
                grid_list = grids_occupied_by_line(f, plg[j][0], plg[i][0])
        grid_lists.append(grid_list)
    return grids_encircled_by_line_grids(grid_lists)


if __name__ == '__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    veh = Vehicle(trajectory=np.array([[-1,10.125,10.125,np.pi/4,0]]))
    # r = np.sqrt(veh.width**2 + (veh.length+veh.wheelbase)**2) / 2
    # n = 2*floor(ceil(8*r)/2)+1
    # print(r,n)
    # vehh = Vehicle(trajectory=np.array([[-1, n/8, n/8, 0.75*np.pi, 0]]))
    K = grids_occupied_by_vehicle(veh)
    # print(K.shape)
    ax1.imshow(K,cmap=plt.cm.gray_r, origin="lower", extent=(0,20.25,0,20.25))

    plg = [(5,5),(10,12),(8,8),(3,2)]
    Env = grids_occupied_by_polygon(plg)
    ax2.imshow(Env,cmap=plt.cm.gray_r, origin="lower", extent=(0,20.25,0,20.25))

    # KK = np.array([[1,1,1,1,2,1,1],[1,2,3,1,1,1,1],[3,2,1,1,1,1,1],[1,1,1,1,2,1,1],[1,1,1,1,2,1,1],[1,1,1,1,2,1,1],[1,1,1,1,2,1,1]])
    Dest = cv2.filter2D(Env, -1, K)
    # Dest = cv2.dilate(Env, K, iterations = 1)
    ax3.imshow(Dest,cmap=plt.cm.gray_r, origin="lower", extent=(0,20.25,0,20.25))

    plt.show()

    #print(veh.vertex[0:2])
    #print(grid_list)
    # r = np.sqrt(veh.length**2+veh.width**2)/2
    # n = 2*floor(ceil(8*r)/2)+1 #pixel size:0.25*0.25
    # print(r)
    # print(n)

    # verts = [tuple(veh.vertex[i]) for i in range(6)]
    # verts.append(verts[0])


    # codes = [Path.MOVETO,
    #     Path.LINETO,
    #     Path.LINETO,
    #     Path.LINETO,
    #     Path.LINETO,
    #     Path.LINETO,
    #     Path.CLOSEPOLY,
    #     ]
    # path = Path(verts, codes)

    # fig = plt.figure() #(figsize=(n/120,n/120),dpi=120)
    # ax = fig.add_subplot(111)
    # patch = patches.PathPatch(path,facecolor='green')
    # ax.add_patch(patch)
    # ax.imshow(R,cmap=plt.cm.gray_r, origin="lower", extent=(0,20,0,20))


    # plt.axis('equal')
    # plt.axis('off')
    # #ax.set_xlim(-0.125*n,0.125*n)
    # #ax.set_ylim(-0.125*n,0.125*n)

    # plt.show()
    # plt.savefig('vehiclemap.png',)

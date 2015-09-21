import numpy as np
import numpy.matlib
# from scipy.integrate import quad
import matplotlib.pyplot as plt
import spiral3


def center_line(q0, k_c, s_g, N):
    # q0: (x0, y0, theta0, k0)
    # k_c: curvature function of path length
    # s_g: total path length
    # N: number of discrete points along path
    delta_s = s_g / (N - 1) # delta_s < 0.1
    s_m = list(np.linspace(0., s_g, N))
    k_m = [k_c(s) for s in s_m]
    theta_m = list(np.zeros(N))
    x_m = list(np.zeros(N))
    y_m = list(np.zeros(N))
    theta_m[0] = q0[2]
    x_m[0] = q0[0]
    y_m[0] = q0[1]
    for i in range(1, N):
        theta_m[i] = theta_m[i - 1] + (k_m[i - 1] + k_m[i]) / 2 * delta_s
    cos_t = [np.cos(t) for t in theta_m]
    sin_t = [np.sin(t) for t in theta_m]
    for i in range(1, N):
        x_m[i] = x_m[i - 1] + (cos_t[i - 1] + cos_t[i]) / 2 * delta_s
        y_m[i] = y_m[i - 1] + (sin_t[i - 1] + sin_t[i]) / 2 * delta_s
    return x_m, y_m, theta_m, k_m


def longitudinal_line(x_m, y_m, theta_m, d):
    #k_d = [(k**-1 - d)**-1 for k in k_m]
    #theta_d = theta_m
    x_d = list(np.zeros(len(x_m)))
    y_d = list(np.zeros(len(x_m)))
    for i in range(len(x_m)):
        x_d[i] = x_m[i] + d*np.cos(theta_m[i]+np.pi/2)
        y_d[i] = y_m[i] + d*np.sin(theta_m[i]+np.pi/2)
    return x_d, y_d # , theta_d, k_d theta_d and k_d may not be needed!


def s_l_coordinate(q_m, s_g, s, l):
    # q_m: (x_m, y_m, theta_m, k_m)
    # s_g: path length
    if 0<=s<=s_g:
        N = len(q_m[0])
        n = int(np.ceil(s/s_g*N) - 1)
        d_s = s_g/N
        d_1 = s - n*d_s
        l1 = d_1 /d_s
        l2 = 1- l1
        k_s = q_m[3][n]*l2 + q_m[3][n+1]*l1
        theta_s = q_m[2][n]*l2 + q_m[2][n+1]*l1
        x_s = q_m[0][n]*l2 + q_m[0][n+1]*l1
        y_s = q_m[1][n]*l2 + q_m[1][n+1]*l1
        ###
        k_sl = (k_s**-1 - l)**-1
        theta_sl = theta_s
        x_sl = x_s + l*np.cos(theta_s+np.pi/2)
        y_sl = y_s + l*np.sin(theta_s+np.pi/2)
        return x_sl, y_sl, theta_sl, k_sl
    else:
        return None


def lateral_line(x_m, y_m, theta_m, s_g, w, s):
    # s_g: path length
    # w: road width
    if 0<=s<=s_g:
        N = len(x_m)
        n = int(np.ceil(s/s_g*N) - 1)
        d_s = s_g/N
        d_1 = s - n*d_s
        l1 = d_1 /d_s
        l2 = 1- l1
        theta_s = theta_m[n]*l2 + theta_m[n+1]*l1
        x_s = x_m[n]*l2 + x_m[n+1]*l1
        y_s = y_m[n]*l2 + y_m[n+1]*l1
        ###
        lateral_bias = list(np.linspace(-w,w,21))
        x = [ x_s + l*np.cos(theta_s+np.pi/2) for l in lateral_bias ]
        y = [ y_s + l*np.sin(theta_s+np.pi/2) for l in lateral_bias ]
        return x, y
    else:
        return None


def grid_width(lane_width, ref_width):
    # ref_width = 0.6
    return lane_width / ( 2 * np.ceil(np.ceil(lane_width / ref_width) / 2) )


if __name__ == '__main__':
    #s_g = 109.61234579137816
    s_g = 211.42520307451244
    #xm, ym, thetam, km = center_line((0,0,0,0.01), lambda s: 0.01-0.000242811907598*s+6.42266190994e-6*s**2-5.3571326595e-8*s**3, s_g, 1001)
    # q0=(0,0,0,0.01), q1=(200,60,pi/6.-0.005)
    xm, ym, thetam, km = center_line((0,0,0,0.01), lambda s: 0.01-0.00038611167606780294*s+4.4656981495145228e-6*s**2-1.4071316854528154e-8*s**3, s_g, 2201)
    sl2xy = lambda s, l:s_l_coordinate((xm,ym,thetam,km),s_g, s, l)
    plt.plot(xm,ym, color='green', linestyle='--', linewidth=2.)
    # print(xm[-1], ym[-1], thetam[-1])
    xu, yu = longitudinal_line(xm, ym, thetam, 5.)
    xl, yl = longitudinal_line(xm, ym, thetam, -5.)
    plt.plot(xu,yu, color='black', linewidth=1.)
    plt.plot(xl,yl, color='black', linewidth=1.)
    ####
    lateral_bias = [-3.75,-2.5,-1.25,0.,1.25,2.5,3.75]
    longitudinal_bias = list(np.arange(0, s_g, 3.75))
    for d in lateral_bias:
        x, y = longitudinal_line(xm, ym, thetam, d)
        plt.plot(x,y, color='black', linewidth=0.3)
    for s in longitudinal_bias:
        x, y = lateral_line(xm, ym, thetam, s_g, 5, s)
        plt.plot(x,y, color='black', linewidth=0.3)
    # 
    '''
    q_s = sl2xy(22.5,2.5)
    plt.plot(q_s[0], q_s[1], '.')
    for d in lateral_bias:
        q_g = sl2xy(45.,d)
        plt.plot(q_g[0], q_g[1], '.')
        p = spiral3.calc_path(q_s, q_g)
        x, y, t, k = center_line(q_s, lambda s: spiral3.k(s,p), p[-1],251)
        plt.plot(x, y, linewidth=0.5)
        #
        q_g2 = sl2xy(90., -2.5)
        plt.plot(q_g2[0], q_g2[1], '.')
        #
        for d in lateral_bias:
            q_g1 = sl2xy(67.5,d)
            plt.plot(q_g1[0], q_g1[1], '.')
            p1 = spiral3.calc_path(q_g, q_g1)
            x, y, t, k = center_line(q_g, lambda s: spiral3.k(s,p1), p1[-1], 251)
            plt.plot(x, y, linewidth=0.5)
            #
            p2 = spiral3.calc_path(q_g1, q_g2)
            x, y, t, k = center_line(q_g1, lambda s: spiral3.k(s,p2), p2[-1], 251)
            plt.plot(x, y, linewidth=0.5)
    '''
    plt.axis('equal')
    plt.savefig('road5.png',dpi=500)
    plt.show()

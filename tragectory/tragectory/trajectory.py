from road import *
import numpy as np
import numpy.matlib
from scipy.integrate import quad
import matplotlib.pyplot as plt
import spiral3

if __name__ == '__main__':
    s_g = 109.61234579137816
    xm, ym, thetam, km = center_line((0,0,0,0.01), lambda s: 0.01-0.000242811907598*s+6.42266190994e-6*s**2-5.3571326595e-8*s**3, s_g, 1001)
    sl2xy = lambda s, l:s_l_coordinate((xm,ym,thetam,km),s_g, s, l)
    plt.plot(xm,ym, color='green', linestyle='--', linewidth=2.)
    # print(x[-1], y[-1], theta[-1])
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
    plt.axis('equal')
    plt.savefig('road4.png',dpi=500)
    plt.show()

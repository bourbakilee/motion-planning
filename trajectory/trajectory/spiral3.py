import numpy as np
import numpy.matlib
from scipy.integrate import quad
import matplotlib.pyplot as plt


def a(p):
    # p = (p0,p1,p2,p3,sg)
    # if p[4] > 0:
    return p[0]
    # else:
    #     return None


def b(p):
    # p = (p0,p1,p2,p3,sg)
    # if p[4] > 0:
    return -(11*p[0]-18*p[1]+9*p[2]-2*p[3])/(2*p[4])
    # else:
        # return None


def c(p):
    # p = (p0,p1,p2,p3,sg)
    # if p[4] > 0:
    return 9*(2*p[0]-5*p[1]+4*p[2]-p[3])/(2*p[4]**2)
    # else:
    #     return None


def d(p):
    # p = (p0,p1,p2,p3,sg)
    # if p[4] > 0:
    return -9*(p[0]-3*p[1]+3*p[2]-p[3])/(2*p[4]**3)
    # else:
    #     return None


def k(s, p):
    # if 0. <= s and s <= p[4]:
    return a(p) + b(p)*s + c(p)*s**2 + d(p)*s**3
    # else:
    #     return None


def dk_ds(s, p):
    # if 0. <= s and s <= p[4]:
    return b(p) + c(p)*s*2 + d(p)*s**2*3
    # else:
    #     return None


def theta(s, p):
    # if 0. <= s and s <= p[4]:
    return a(p)*s + b(p)*s**2/2 + c(p)*s**3/3 + d(p)*s**4/4
    # else:
    #     return None


def x(s, p):
    # if 0. <= s and s <= p[4]:
    return quad(lambda t: np.cos(theta(t, p)), 0, s)[0]
    # else:
    #     return None


def y(s, p):
    # if 0. <= s and s <= p[4]:
    return quad(lambda t: np.sin(theta(t, p)), 0, s)[0]
    # else:
    #     return None


def s(t, p):
    # p : (v0, a0, q2, tf)
    # if 0. <= t and t <= p[3]:
    return p[0]*t + p[1]*t**2/2 + p[2]*t**3/3
    # else:
    #     return None


def v(t, p):
    # p : (v0, a0, q2, tf)
    if 0. <= t and t <= p[3]:
        return p[0] + p[1]*t + p[2]*t**2
    else:
        return None


def dk_dt(t, p1, p2):
    # p1 : parameters of path
    # p2 : parameters of velocity
    return v(t, p2)*dk_ds(s(t, p2), p1)
##################################


def q_g(bd_con):
    # bd_con : boundary condition, (k0, x1, y1, theta1, k1)
    # pp : matrix[p1; p2; sg]
    # return : matrix[x;y;theta]
    return lambda pp: np.matrix([[x(pp[2, 0], [bd_con[0], pp[0, 0], pp[1, 0], bd_con[4], pp[2, 0]])],
        [y(pp[2, 0], [bd_con[0], pp[0, 0], pp[1, 0], bd_con[4], pp[2, 0]])],
        [theta(pp[2, 0], [bd_con[0], pp[0, 0], pp[1, 0], bd_con[4], pp[2, 0]])]])


def Jac(pp, bd_con):
    # pp : matrix[p1; p2; sg]
    # bd_con : boundary condition, (k0, x1, y1, theta1, k1)
    # return : matrix[dx/dp1,dx/dp2,dx/dsg;dy...;dtheta...]
    q = q_g(bd_con)
    J = np.matlib.zeros((3, 3)) # Matrix
    J[:, 0] = ( q(np.matrix([[pp[0, 0] + 0.01], [pp[1, 0]], [pp[2, 0]]])) - q(np.matrix([[pp[0, 0] - 0.01], [pp[1, 0]], [pp[2, 0]]])) )/0.02
    J[:, 1] = ( q(np.matrix([[pp[0, 0]], [pp[1, 0] + 0.01], [pp[2, 0]]])) - q(np.matrix([[pp[0, 0]], [pp[1, 0] - 0.01], [pp[2, 0]]])) )/0.02
    J[:, 2] = ( q(np.matrix([[pp[0, 0]], [pp[1, 0]], [pp[2, 0] + 0.1]])) - q(np.matrix([[pp[0, 0]], [pp[1, 0]], [pp[2, 0] - 0.1]])) )/0.2
    return J


def delta_q(pp, bd_con):
    # pp : matrix[p1; p2; sg]
    # bd_con : boundary condition, (k0, x1, y1, theta1, k1)
    # return : matrix[delta_x; delta_y; delta_theta]
    q1 = np.matrix([[bd_con[1]], [bd_con[2]], [bd_con[3]]])
    q = q_g(bd_con)
    q_sg = q(pp)
    return q1 - q_sg


def dist(q0, q1):
    # q = matrix[x; y; theta]
    k_max = 0.2
    dtheta = np.mod(abs(q1[2, 0] - q0[2, 0]), 2*np.pi)
    return np.sqrt((q1[0, 0] - q0[0, 0])**2 + (q1[1, 0] - q0[1, 0])**2) + min(dtheta, 2*np.pi - dtheta) / k_max


def opt_path(bd_con, init_val):
    # bd_con : boundary condition, (k0, x1, y1, theta1, k1)
    # init_val: np.matrix - 3X1
    # return: np.matrix - 3X1
    norm = lambda q:dist(np.matlib.zeros((3,1)), q)
    pp = init_val
    # pp = np.matrix([
        # [(2.*bd_con[0]+bd_con[4])/3.],
        # [(bd_con[0]+2.*bd_con[4])/3.],
        # [norm(np.matrix([[bd_con[1]], [bd_con[2]], [bd_con[3]]]))]
        # ]) # initial guess value, can be replaced by lookup-table
    norm_dq = 1.
    eps = 1.e-8
    # N = 100
    while norm_dq > eps:
        J = Jac(pp, bd_con)
        dq = delta_q(pp, bd_con)
        pp += J**-1*dq
        norm_dq = norm(dq)
        # print(pp)
        # print(norm_q)
    return pp


def calc_path(q0, q1, init_val=None):
    # q : (x, y, theta, k)
    # return : (p0,p1,p2,p3,sg)
    cc = np.cos(q0[2])
    ss = np.sin(q0[2])
    x_r = (q1[0]-q0[0])*cc + (q1[1]-q0[1])*ss
    y_r = -(q1[0]-q0[0])*ss + (q1[1]-q0[1])*cc
    theta_r = np.mod(q1[2]-q0[2], 2*np.pi)
    bd_con = (q0[3], x_r, y_r, theta_r, q1[3])
    norm = lambda q:dist(np.matlib.zeros((3,1)), q)
    if init_val is None:
        init_val =  np.matrix([
        [(2.*bd_con[0]+bd_con[4])/3.],
        [(bd_con[0]+2.*bd_con[4])/3.],
        [norm(np.matrix([[bd_con[1]], [bd_con[2]], [bd_con[3]]]))]
        ]) # initial guess value, can be replaced by lookup-table
    pp = opt_path(bd_con, init_val)
    p = (q0[3], pp[0, 0], pp[1, 0], q1[3], pp[2, 0])
    return p


def calc_velocity(v0, a0, vf, sf):
    # a0 selection: a0 = amax*(1-v0/vmax)
    eps = 1.e-6
    if abs(a0) > eps:
        delta = (2*v0+vf)**2 + 6*a0*sf
        if delta >= 0:
            tf = (np.sqrt(delta)-2*v0-vf)/a0
            if tf > eps:
                q2 = (vf-v0)/tf**2 - a0/tf
                return (v0, a0, q2, tf)
            else:
                return None
        else:
            return None
    elif abs(a0) <= eps and (v0 > eps or vf > eps):
        tf = 3*sf/(2*v0+vf)
        q2 = (vf-v0)/tf**2
        return (v0, a0, q2, tf)
    else:
        return None


def eval_path(p, plims):
    # p : parameters of path (p0,p1,p2,p3,sg)
    # plims : (k_m, v_max, a_max, a_min, dk_m)
    path = True
    bb, cc, dd = b(p), c(p), d(p)
    delta = cc**2+3*bb*dd
    if delta >= 0:
        s1 = (-cc + np.sqrt(delta))/(3*dd)
        s2 = (-cc - np.sqrt(delta))/(3*dd)
        if ((0 <= s1 <= p[4] and abs(k(s1, p)) <= plims[0]) or
                (0 <= s2 <= p[4] and abs(k(s2, p)) <= plims[0])):
            path = False
    return path


def eval_velocity(p1, p2, plims):
    # p1 : parameters of path (p0,p1,p2,p3,sg)
    # p2 : parameters of velocity (q0,q1,q2,tg)
    # plims : (k_m, v_max, a_max, a_min, dk_m)
    eps = 1.e-6
    velocity = True
    a_g = p2[1]+2*p2[2]*p2[3]
    if a_g > plims[2] or a_g < plims[3]:
        velocity = False
    elif p2[2] > eps:
        tm = -p2[1]/(2*p2[2])
        if 0 <= tm <= p2[3] and v(tm, p2) > plims[1]:
            velocity = False
    elif p2[2] < -eps:
        tm = -p2[1]/(2*p2[2])
        if 0 <= tm <= p2[3] and v(tm, p2) < 0.:
            velocity = False
    else:
        def f(t): return dk_dt(t, p1, p2)
        ttt = np.linspace(0., p2[3], np.floor(10*tt) + 2)
        vv = [t for t in ttt if abs(f(t)) > plims[4]]
        if vv:
            velocity = False
    return velocity


if __name__ == '__main__':
    # q0 = [1., 1., np.pi/6, 0.01]
    # q1 = [100.,50.,np.pi/3,-0.01]
    # p1 = calc_path(q0,q1)
    # plims = (0.2, 12., 1., -5., 0.1)
    # v0,a0,vf = 6., 0.5, 8.
    # p2 = calc_velocity(v0,a0,vf,p1[4])
    # print(p1)
    # print(eval_path(p1,plims))
    # print(p2)
    # print(eval_velocity(p1,p2,plims))

    bd_con = (0.01, 100, 40, np.pi/6., -0.01)
    q0 = (0, 0, 0, bd_con[0])
    q1 = (bd_con[1], bd_con[2], bd_con[3], bd_con[4])
    p = calc_path(q0,q1)
    #p = [bd_con[0], pp[0, 0], pp[1, 0], bd_con[4], pp[2, 0]]
    sp = np.linspace(0, p[4], 101)
    xx = [x(ss, p) for ss in sp]
    yy = [y(ss, p) for ss in sp]
    tt = [theta(ss, p) for ss in sp]
    print(p)
    print(a(p),b(p),c(p),d(p))
    print(xx[-1])
    print(yy[-1])
    print(tt[-1])
    plt.plot(xx, yy)
    plt.axis('equal')
    plt.show()

import numpy as np
import numpy.matlib
from scipy.integrate import quad
import matplotlib.pyplot as plt


def a(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return p[0]
    else:
        return None


def b(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return p[1]
    else:
        return None


def c(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return p[2]/2
    else:
        return None


def d(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return -(575*p[0]-648*p[3]+81*p[4]-8*p[5]+170*p[1]*p[6]+22*p[2]*p[6]**2)/(8*p[6]**3)
    else:
        return None


def e(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return 9*(37*p[0]-45*p[3]+9*p[4]-p[5]+10*p[1]*p[6]+p[2]*p[6]**2)/(2*p[6]**4)
    else:
        return None


def f(p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if p[6] > 0:
        return -9*(85*p[0]-108*p[3]+27*p[4]-4*p[5]+22*p[1]*p[6]+2*p[2]*p[6]**2)/(8*p[6]**5)
    else:
        return 


def k(s, p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if 0. <= s and s <= p[6]:
        return a(p) + b(p)*s + c(p)*s**2 + d(p)*s**3 + e(p)*s**4 + f(p)*s**5
    else:
        return None


def dk_ds(s, p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if 0. <= s and s <= p[6]:
        return b(p) + c(p)*s*2 + d(p)*s**2*3 + e(p)*s**3*4 + f(p)*s**4*5
    else:
        return None


def theta(s, p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if 0. <= s and s <= p[6]:
        return a(p)*s + b(p)*s**2/2 + c(p)*s**3/3 + d(p)*s**4/4 + e(p)*s**5/5 + f(p)*s**6/6
    else:
        return None


def x(s, p):
    # p = (p0,p1,p2,p3,p4,p5,sg)
    if 0. <= s and s <= p[6]:
        return quad(lambda t: np.cos(theta(t, p)), 0, s)[0]
    else:
        return None


def y(s, p):
    if 0. <= s and s <= p[6]:
        return quad(lambda t: np.sin(theta(t, p)), 0, s)[0]
    else:
        return None


def q_g(bd_con):
    # bd_con : boundary condition, (k0, dk0, ddk0, x1, y1, theta1, k1)
    # pp : matrix[p3; p4; sg]
    # return : matrix[x;y;theta]
    return lambda pp: np.matrix([[x(pp[2, 0], ([bd_con[0], bd_con[1], bd_con[2], pp[0, 0], pp[1, 0], bd_con[6], pp[2, 0]]))],
        [y(pp[2, 0], ([bd_con[0], bd_con[1], bd_con[2], pp[0, 0], pp[1, 0], bd_con[6], pp[2, 0]]))],
        [theta(pp[2, 0], ([bd_con[0], bd_con[1], bd_con[2], pp[0, 0], pp[1, 0], bd_con[6], pp[2, 0]]))]])


def Jac(pp, bd_con):
    # pp : matrix[p3; p4; sg]
    # bd_con : boundary condition, (k0, dk0, ddk0, x1, y1, theta1, k1)
    # return : matrix[dx/dp3,dx/dp4,dx/dsg;dy...;dtheta...]
    q = q_g(bd_con)
    J = np.matlib.zeros((3, 3)) # Matrix
    J[:, 0] = ( q(np.matrix([[pp[0, 0] + 0.005], [pp[1, 0]], [pp[2, 0]]])) - q(np.matrix([[pp[0, 0] - 0.005], [pp[1, 0]], [pp[2, 0]]])) )/0.01
    J[:, 1] = ( q(np.matrix([[pp[0, 0]], [pp[1, 0] + 0.005], [pp[2, 0]]])) - q(np.matrix([[pp[0, 0]], [pp[1, 0] - 0.005], [pp[2, 0]]])) )/0.01
    J[:, 2] = ( q(np.matrix([[pp[0, 0]], [pp[1, 0]], [pp[2, 0] + 0.1]])) - q(np.matrix([[pp[0, 0]], [pp[1, 0]], [pp[2, 0] - 0.1]])) )/0.2
    return J


def delta_q(pp, bd_con):
    # pp : matrix[p3; p4; sg]
    # bd_con : boundary condition, (k0, dk0, ddk0, x1, y1, theta1, k1)
    # return : matrix[delta_x; delta_y; delta_theta]
    q1 = np.matrix([[bd_con[3]], [bd_con[4]], [bd_con[5]]])
    q = q_g(bd_con)
    q_sg = q(pp)
    return q1 - q_sg

def dist(q0, q1):
    # q = matrix[x; y; theta]
    k_max = 0.2
    dtheta = np.mod(abs(q1[2, 0] - q0[2, 0]), 2*np.pi)
    return np.sqrt((q1[0, 0] - q0[0, 0])**2 + (q1[1, 0] - q0[1, 0])**2) + min(dtheta, 2*np.pi - dtheta) / k_max


def opt_path(bd_con):
    # bd_con : boundary condition, (k0, dk0, ddk0, x1, y1, theta1, k1)
    # pp : matrix[p3; p4; sg]
    norm = lambda q:dist(np.matlib.zeros((3,1)), q)
    pp = np.matrix([
        [(2.*bd_con[0]+bd_con[6])/3.],
        [(bd_con[0]+2.*bd_con[6])/3.],
        [norm(np.matrix([[bd_con[3]], [bd_con[4]], [bd_con[5]]]))]
        ]) # initial guess value, can be replaced by lookup-table
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


def calc_path(q0, q1, dk0, ddk0):
    # q : (x, y, theta, k)
    # return : (p0,p1,p2,p3,p4,p5,sg)
    cc = np.cos(q0[2])
    ss = np.sin(q0[2])
    x_r = (q1[0]-q0[0])*cc + (q1[1]-q0[1])*ss
    y_r = -(q1[0]-q0[0])*ss + (q1[1]-q0[1])*cc
    theta_r = np.mod(q1[2]-q0[2], 2*np.pi)
    bd_con = (q0[3], dk0, ddk0, x_r, y_r, theta_r, q1[3])
    pp = opt_path(bd_con)
    p = (q0[3], dk0, ddk0, pp[0, 0], pp[1, 0], q1[3], pp[2, 0])
    return p


if __name__ == '__main__':
    q0 = [1., 1., np.pi/6, 0.01]
    q1 = [100.,50.,np.pi/3,-0.01]
    p1 = calc_path(q0,q1, 0., 0.)
    # print(p)
    # print(a(p), b(p), c(p), d(p), e(p), f(p))
    import road
    xx1 ,yy1, tt1, kk1 = road.center_line(q0, lambda s:k(s,p1), p1[-1], 2001)
    import spiral3
    p2 = spiral3.calc_path(q0, q1)
    xx2 ,yy2, tt2, kk2 = road.center_line(q0, lambda s:spiral3.k(s,p2), p2[-1], 2001)
    plt.plot(xx1,yy1, label='spiral5')
    plt.plot(xx2, yy2, label='spiral3')
    plt.legend()
    plt.axis('equal')
    plt.show()

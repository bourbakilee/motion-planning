import numpy as np

def spiral_calc(fun, length,q=None,ref_delta_s=0.1):
    # fun: k(s)
    # q0=(0,0,0)
    N = int(np.ceil(length/ref_delta_s))
    line = np.zeros((N,4))
    delta_s = length / (N-1)
    s_list = np.linspace(0,length,N)
    line[:,3] = fun(s_list) # k
    d_theta = (line[0:N-1,3] + line[1:N,3])/2 * delta_s
    for i in range(1,N):
        line[i,2] = line[i-1,2] + d_theta[i-1] # theta
    cos_t = np.cos(line[:,2])
    sin_t = np.sin(line[:,2])
    d_x = (cos_t[0:N-1] + cos_t[1:N])/2 * delta_s
    d_y = (sin_t[0:N-1] + sin_t[1:N])/2 * delta_s
    for i in range(1,N):
        line[i,0] = line[i-1,0] + d_x[i-1] # x
        line[i,1] = line[i-1,1] + d_y[i-1] # y
    if q is not None:
        sin_x = np.sin(q[2])*line[:,0]
        cos_x = np.cos(q[2])*line[:,0]
        sin_y = np.sin(q[2])*line[:,1]
        cos_y = np.cos(q[2])*line[:,1]
        line[:,0] = q[0] + cos_x - sin_y
        line[:,1] = q[1] + sin_x + cos_y
        line[:,2] = np.mod((line[:,2] + q[2]), 2*np.pi)
    return line
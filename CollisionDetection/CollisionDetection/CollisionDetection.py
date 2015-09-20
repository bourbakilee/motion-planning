import numpy as np
from scipy.fftpack import fft2, ifft2
import pylab as pl

N=150 # map size
o = np.zeros((N,N)) # obstacles
o[0, :] = 1
o[N-1, :] = 1
o[:, 0] = 1
o[:, N-1] = 1
o[60:90, 60:90] = 1
r = np.zeros((N, N)) # robot
r[0:5, 0:5] = 1
r[N-4:N, 0:5]= 1
r[0:5, N-4:N] = 1
r[N-4:N, N-4:N] = 1
O = fft2(o)
R = fft2(r)
C = O*R # convolution
c = ifft2(C)
pl.imshow(np.real(c), cmap=pl.cm.gray_r)

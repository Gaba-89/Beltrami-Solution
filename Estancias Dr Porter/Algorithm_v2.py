import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
import copy
import time
#import multiprocessing as mp



# Update a triangle.
def newW1(nup,z1,z2,z3,w1,w2,w3):
    p=((z1-z3)/(z1-z2)+nup*(z1.conjugate()-z3.conjugate())/(z1.conjugate()-z2.conjugate()))/(1+nup)
    p0, p1 = p.real, p.imag
    u1, v1 =  w1.real, w1.imag
    u2, v2 =  w2.real, w2.imag
    u3, v3 =  w3.real, w3.imag
    abspsq = p0*p0 + p1*p1
    delta = 2*(abspsq-p0+1)
    a0 = ( -(u1-2*u2+u3)-p0*(u1+u2-2*u3)-p1*(v1+v2-2*v3) )/delta
    a1 = ( -(v1-2*v2+v3)-p0*(v1+v2-2*v3)+p1*(u1+u2-2*u3) )/delta
    b0 = (  (u1+u3)-p0*(u2+u3)+p1*(v2-v3)+abspsq*(u1+u2) )/delta
    b1 = (  (v1+v3)-p0*(v2+v3)-p1*(u2-u3)+abspsq*(v1+v2) )/delta
    b = b0 + b1*1j
    a = a0 + a1*1j

    return b, a + b, a*p+b

# for drawing fig. This should be parallelized.
def drawmesh(w,N,M):
    plt.figure(figsize=(10*M/(M+N), 10*N/(M+N)), dpi=120)
    for j in jlist:
        plt.plot(w[j].real,w[j].imag,'b-', lw=0.5)
    for k in klist:
        wtemp = np.zeros((len(jlist)),dtype=complex)
        wtemp = np.array([w[j][k] for j in jlist])
        ccode = 'r-' if (k%2)==0 else 'g-'
        plt.plot(wtemp.real,wtemp.imag, ccode, lw=0.5)
        if k>0:
            wtemp = np.array([w[j][k-(j%2)] for j in jlist])
            plt.plot(wtemp.real,wtemp.imag, 'k-', lw=0.5)
    plt.show()



def first_phase(z,w,N,M):
    sigma = np.zeros(((N-1)*M*2,3),dtype=complex)
    n = 0
    l = (N-1)*2
    for j in range(M):
      for k in range(N-1):
        if j % 2 == 0:
          # Beltrami Coefficient of z-dom 
          def mu(z):
            return 0.0

# Beltrami Coefficient on Z-dom.
          def nu(z):
            return mu(cmath.exp(z))*cmath.exp(-2*1j*z.imag)
          
          w1t, w2t, w3t = newW1(nu(z[-1-j][k]),
                               z[-j][k],z[-j][k+1],z[-1-j][k],
                               w[-j][k],w[-j][k+1],w[-1-j][k])
          sigma[n][0] = w1t
          sigma[n][1] = w2t
          sigma[n][2] = w3t

          w1b, w2b, w3b = newW1(nu(z[-j][k+1]),
                                z[-j-1][k],z[-j-1][k+1],z[-j][k+1],
                               w[-j-1][k],w[-j-1][k+1],w[-j][k+1])
          sigma[n+(N-1)][0] = w1b
          sigma[n+(N-1)][1] = w2b
          sigma[n+(N-1)][2] = w3b
          n += 1
        else:
          
                    # Beltrami Coefficient of z-dom 
          value = np.linspace(0.0,0.9,4)
          print(value[k])
          def mu(z):
            return value[k]

# Beltrami Coefficient on Z-dom.
          def nu(z):
            return mu(cmath.exp(z))*cmath.exp(-2*1j*z.imag)
          w1t, w2t, w3t = newW1(nu(z[-j][k]),
                                z[-1-j][k],z[-1-j][k+1],z[-j][k],
                                w[-1-j][k],w[-1-j][k+1],w[-j][k])
          sigma[l+(N-1)][0] = w1t
          sigma[l+(N-1)][1] = w2t
          sigma[l+(N-1)][2] = w3t

          w1b, w2b, w3b = newW1(nu(z[-j-1][k+1]),
                                z[-j][k],z[-j][k+1],z[-j-1][k+1],
                               w[-j][k],w[-j][k+1],w[-j-1][k+1])
          sigma[l][0] = w1b
          sigma[l][1] = w2b
          sigma[l][2] = w3b
          l += 1
      if j % 2 == 0:
        n += (N-1)*3
      else:
        l += (N-1)*3

    return sigma




def plot_triangles(triangles):
  edge_colors = ["black","red","blue"]
  for verts in triangles:
      x = [z.real for z in verts] + [verts[0].real]
      y = [z.imag for z in verts] + [verts[0].imag]
      for i in range(3):
        x_edge = [x[i], x[(i+1)%3]]
        y_edge = [y[i], y[(i+1)%3]]
        plt.plot(x_edge, y_edge, color=edge_colors[i%len(edge_colors)], linewidth=.9)
# Input
N,M =5,5

# Initialize Z-Vertices and W-Vertices
klist = range(N)
jlist = range(-M,1)
r = math.sqrt(3.)*math.pi/N
z = np.zeros((len(jlist),len(klist)),dtype=complex)
w = np.zeros((len(jlist),len(klist)),dtype=complex)
for j in jlist:
    for k in klist:
        z[j][k] = r*j + 2*math.pi*1j*(k+(j%2)/2.)/N
w = copy.deepcopy(z)

drawmesh(z,N,M)
mesh_sigma = first_phase(z,w,N,M)

plot_triangles(mesh_sigma)
plt.gca().set_aspect("equal", adjustable="box")
plt.show()
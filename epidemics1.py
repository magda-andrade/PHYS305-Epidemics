#Magda Andrade
#Epidemic Program
#April 30


import matplotlib.pyplot as plt
import numpy as np



def RK4(t, y, h, sirRHS, *P):

    thalf = t + 0.5*h
    k1 = h * sirRHS(t, y, *P)
    k2 = h * sirRHS(thalf, y + 0.5*k1, *P)
    k3 = h * sirRHS(thalf, y + 0.5*k2, *P)
    k4 = h * sirRHS(t + h, y + k3, *P)
    return y + (k1 + 2*k2 + 2*k3 + k4)/6

def sirRHS(t, y, *P):   
    beta, gamma = P
    S = y[0]
    y[0]=N -y[1] -y[2]
    I = y[1]
    R = y[2]
    S,I,R = y
    dSdt = -beta * I * S / N
    dIdt = (beta * I * S / N) - (gamma * I)
    dRdt = gamma * I
    
    return np.array([ dSdt, dIdt, dRdt ])

def odeSolve(t0, y0, tmax, h, sirRHS, method, *P):
   
    t = np.arange(t0,tmax+h,h)
    ntimes,  = t.shape

    if type(y0) in [int, float]: 
        neqn = 1
        y = np.zeros( ntimes )
        #neqn, = y0.shape
        #y = np.zeros( (ntimes, neqn) )
    else:                         
        neqn, = y0.shape
        y = np.zeros( (ntimes, neqn) )
 
    y[0] = y0

    for i in range(0,ntimes-1):
        y[i+1] = method(t[i], y[i], h, sirRHS, *P)

    return t,y

beta = 0.5
gamma = 0.3
array = (beta,gamma)
N=100
h = 0.5
t0 = 0.0
y0 = np.array([ 0.6 , 0.3, 0.1 ])
tmax = 150


t, y = odeSolve(t0, y0, tmax, h, sirRHS, RK4, *array)

fig,ax = plt.subplots()
plt.title('When R0 = Rcrit', fontdict=None, loc='center',pad=None )
ax.plot(t,y[:,0],'b', label='Susceptible')
ax.plot(t,y[:,1],'r', label='Infected')
ax.plot(t,y[:,2],'c', label='recovered')
ax.set_xlabel('time')
ax.set_ylabel('population')
ax.legend()

plt.tight_layout()
plt.show()

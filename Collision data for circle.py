import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def collision_data_nospin(N, x0, y0, alpha0):
    
    t = np.zeros(N) # times at collision
    x = np.zeros(N) # x-values at collision
    y = np.zeros(N) # y-values at collision
    phi = np.zeros(N) # phi values at collision
    alpha = np.zeros(N) # alpha values at collision
    b = np.zeros(N) # b values at collision

    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    phi[0] = 0
    alpha[0] = alpha0
    b[0] = x[0]*np.cos(alpha[0]) + y[0]*np.sin(alpha[0])
    c = x[0]**2 + y[0]**2 - 1 #c is only needed for one time-step
    
    # Update formulae for first collision
    t[1] = -b[0] + np.sqrt(b[0]**2 - c)
    x[1] = x[0] + t[1]*np.cos(alpha[0])
    y[1] = y[0] + t[1]*np.sin(alpha[0])
    phi[1] = np.arctan2(y[1], x[1])
    alpha[1] = 2*phi[1] - alpha[0] + np.pi
    b[1] = x[1]*np.cos(alpha[1]) + y[1]*np.sin(alpha[1])
    
    # Update formulae for remaining collisions
    for i in range(2, N):
        t[i] = -2*b[i - 1]
        x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
        y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
        phi[i] = np.arctan2(y[i], x[i])
        alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        b[i] = x[i]*np.cos(alpha[i]) + y[i]*np.sin(alpha[i])
    
    return x, y, alpha

def collision_data(N, x0, y0, alpha0, u0, MI_coeff):
    
    vx0 = np.cos(alpha0)
    vy0 = np.sin(alpha0)
    
    t = np.zeros(N) # times at collision
    x = np.zeros(N) # x-values at collision
    y = np.zeros(N) # y-values at collision
    vx = np.zeros(N) # x-component of veolcity at collision
    vy = np.zeros(N) # y-component of velocity at collision
    u = np.zeros(N) # spin values at collision
    b = np.zeros(N) # b values at collision
    c = np.zeros(N) # c values at collision
    
    t[0] = 0
    x[0] = x0
    y[0] = y0
    vx[0] = vx0
    vy[0] = vy0
    u[0] = u0
    b[0] = 0
    c[0] = 0
    
    for i in range(1,N):
        b[i] = (x[i - 1]*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
        c[i] = (x[i - 1]**2 + y[i - 1]**2 - 1)/(vx[i - 1]**2 + vy[i - 1]**2)
        t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
        x[i] = x[i - 1] + t[i]*vx[i - 1]
        y[i] = y[i - 1] + t[i]*vy[i - 1]
        
        T = np.array([-y[i],x[i]])
        n = np.array([-x[i],-y[i]])
        v = np.array([vx[i - 1], vy[i - 1]])
        
        vT = np.dot(v, T)
        vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
        vperp = -np.dot(v,n)
        
        vnew = vparr*T + vperp*n
        
        u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
        vx[i] = vnew[0]
        vy[i] = vnew[1]
        
    return x, y, u, vx, vy

x_spin, y_spin, u_spin, vx_spin, vy_spin = collision_data(50,0.5,-0.1,np.pi/4,1,1/2)

#x1, y1, alpha1 = collision_data_nospin(50, 0.5, -0.1, np.pi/4)

# May be uncommented to save collision data:
#np.savetxt('circle_edges.txt',np.transpose(np.array([x_spin,y_spin,vx_spin,vy_spin,u_spin])))
  
# Plotting:
thetavals = np.linspace(0, 2*np.pi, 100)
'''
plt.scatter(x1,y1,c='blue')
plt.plot(x1,y1,c='green', label='Unspinable')
'''
plt.scatter(x_spin,y_spin,c='purple')
plt.plot(x_spin,y_spin,c='red')

plt.plot(np.cos(thetavals), np.sin(thetavals), color='black')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()
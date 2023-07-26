import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Width of rectangle:
L = 2

def collision_data_nospin(N, x0, y0, alpha0):
    
    t = np.zeros(N) # times at collision
    x = np.zeros(N) # x-values at collision
    y = np.zeros(N) # y-values at collision
    alpha = np.zeros(N) # alpha values at collision
    theta1 = np.zeros(N) # theta1 values at collision
    theta2 = np.zeros(N) # theta2 values at collision
    theta3 = np.zeros(N) # theta3 values at collision
    theta4 = np.zeros(N) # theta4 values at collision
    
    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    alpha[0] = alpha0
    theta1[0] = np.arctan2(1 - y[0], L - x[0])
    theta2[0] = np.arctan2(1 - y[0], -L - x[0])
    theta3[0] = np.arctan2(-1 - y[0], -L - x[0])
    theta4[0] = np.arctan2( -1 - y[0], L - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (1 - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = 1
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            t[i] = (-L - x[i - 1])/np.cos(alpha[i - 1])
            x[i] = -L
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            alpha[i] = np.pi - alpha[i - 1]
        if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
            t[i] = (-1 - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = -1
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            t[i] = (L - x[i - 1])/np.cos(alpha[i - 1])
            x[i] = L
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            alpha[i] = np.pi - alpha[i - 1]
        theta1[i] = np.arctan2(1 - y[i], L - x[i])
        theta2[i] = np.arctan2(1 - y[i], -L - x[i])
        theta3[i] = np.arctan2(-1 - y[i], -L - x[i])
        theta4[i] = np.arctan2( -1 - y[i], L - x[i])
        
    return x, y, alpha

def collision_data(N, x0, y0, alpha0, u0, MI_coeff):
    
    vx0 = np.cos(alpha0)
    vy0 = np.sin(alpha0)
    
    t = np.zeros(N) # times at collision
    x = np.zeros(N) # x-values at collision
    y = np.zeros(N) # y-values at collision
    vx = np.zeros(N) # x-component of veolcity at collision
    vy = np.zeros(N) # y-component of veolcity at collision
    u = np.zeros(N) # spin values at collision
    alpha = np.zeros(N) # alpha values at collision
    theta1 = np.zeros(N) # theta1 values at collision
    theta2 = np.zeros(N) # theta2 values at collision
    theta3 = np.zeros(N) # theta3 values at collision
    theta4 = np.zeros(N) # theta4 values at collision
    
    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    vx[0] = vx0
    vy[0] = vy0
    u[0] = u0
    alpha[0] = alpha0
    theta1[0] = np.arctan2(1 - y[0], L - x[0])
    theta2[0] = np.arctan2(1 - y[0], -L - x[0])
    theta3[0] = np.arctan2(-1 - y[0], -L - x[0])
    theta4[0] = np.arctan2( -1 - y[0], L - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (1 - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = 1
            
            vT = -vx[i - 1]
            vn = -vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = -vparr
            vy[i] = -vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
            
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            t[i] = (-L - x[i - 1])/vx[i - 1]
            x[i] = -L
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            vT = -vy[i - 1]
            vn = vx[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vperp
            vy[i] = -vparr
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
            t[i] = (-1 - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = -1
            
            vT = vx[i - 1]
            vn = vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vparr
            vy[i] = vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
            
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            t[i] = (L - x[i - 1])/vx[i - 1]
            x[i] = L
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            vT = vy[i - 1]
            vn = -vx[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = -vperp
            vy[i] = vparr
            alpha[i] = np.arctan2(vy[i],vx[i])
        theta1[i] = np.arctan2(1 - y[i], L - x[i])
        theta2[i] = np.arctan2(1 - y[i], -L - x[i])
        theta3[i] = np.arctan2(-1 - y[i], -L - x[i])
        theta4[i] = np.arctan2( -1 - y[i], L - x[i])
        
    return x, y, alpha, vx, vy, u

x1, y1, alpha1 = collision_data_nospin(50, 0.5, -0.5, np.pi/4 + 0.33)
x_spin, y_spin, alpha_spin, vx_spin, vy_spin, u_spin = collision_data(50, -0.5*L, 0.25, -np.pi/4 + 0.2, 0, 1/2)

# May be uncommented to save collision data:
#np.savetxt('rectangle_edges.txt',np.transpose(np.array([x_spin,y_spin,vx_spin,vy_spin,u_spin])))
'''
plt.scatter(x1,y1)
plt.plot(x1,y1, c='green')
'''
plt.scatter(x_spin,y_spin, c='purple')
plt.plot(x_spin,y_spin, c='red')
plt.gca().add_patch(Rectangle((-L,-1),2*L,2,
                    edgecolor='black',
                    facecolor='none'))
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()
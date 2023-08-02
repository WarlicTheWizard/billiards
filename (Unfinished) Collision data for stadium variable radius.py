import numpy as np
import matplotlib.pyplot as plt

# Width of rectangular part of the stadium:
L = 1

# Radius of the circular part of the stadium
R = 1

def collision_data_nospin(N, x0, y0, alpha0):
    t = np.zeros(N) # times at collision
    x = np.zeros(N) # x-values at collision
    y = np.zeros(N) # y-values at collision
    phi = np.zeros(N) # phi values at collision
    alpha = np.zeros(N) # alpha values at collision
    b = np.zeros(N) # b values at collision
    c = np.zeros(N) # c values at collision
    theta1 = np.zeros(N) # theta1 values at collision
    theta2 = np.zeros(N) # theta2 values at collision
    theta3 = np.zeros(N) # theta3 values at collision
    theta4 = np.zeros(N) # theta4 values at collision

    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    phi[0] = 0
    alpha[0] = alpha0
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(R - y[0], L - x[0])
    theta2[0] = np.arctan2(R - y[0], -L - x[0])
    theta3[0] = np.arctan2(-R - y[0], -L - x[0])
    theta4[0] = np.arctan2( -R - y[0], L - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (R - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = R
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            b[i] = (x[i - 1] + L)*np.cos(alpha[i - 1]) + y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] + L)**2 + y[i - 1]**2 - R**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] + L)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
            t[i] = (-R - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = -R
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            b[i] = (x[i - 1] - L)*np.cos(alpha[i - 1]) + y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] - L)**2 + y[i - 1]**2 - R**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] - L)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        theta1[i] = np.arctan2(R - y[i], L - x[i])
        theta2[i] = np.arctan2(R - y[i], -L - x[i])
        theta3[i] = np.arctan2(-R - y[i], -L - x[i])
        theta4[i] = np.arctan2( -R - y[i], L - x[i])
        
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
    phi = np.zeros(N) # phi values at collision
    alpha = np.zeros(N) # alpha values at collision
    b = np.zeros(N) # b values at collision
    c = np.zeros(N) # c values at collision
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
    phi[0] = 0
    alpha[0] = alpha0
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(R - y[0], L - x[0])
    theta2[0] = np.arctan2(R - y[0], -L - x[0])
    theta3[0] = np.arctan2(-R - y[0], -L - x[0])
    theta4[0] = np.arctan2( -R - y[0], L - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (R - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = R
            
            vT = -vx[i - 1]
            vn = -vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = -vparr
            vy[i] = -vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] + L)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] + L)**2 + y[i - 1]**2 - R**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] + L])
            n = np.array([-(x[i] + L),-y[i]])
            v = np.array([vx[i - 1], vy[i - 1]])
            
            vT = np.dot(v, T)
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -np.dot(v,n)
            
            vnew = vparr*T + vperp*n
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vnew[0]
            vy[i] = vnew[1]
            alpha[i] = np.arctan2(vy[i], vx[i])
        if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
            t[i] = (-R - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = -R
            
            vT = vx[i - 1]
            vn = vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vparr
            vy[i] = vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] - L)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] - L)**2 + y[i - 1]**2 - R**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] - L])
            n = np.array([-(x[i] - L),-y[i]])
            v = np.array([vx[i - 1], vy[i - 1]])
            
            vT = np.dot(v, T)
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -np.dot(v,n)
            
            vnew = vparr*T + vperp*n
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vnew[0]
            vy[i] = vnew[1]
            alpha[i] = np.arctan2(vy[i], vx[i])
        theta1[i] = np.arctan2(R - y[i], L - x[i])
        theta2[i] = np.arctan2(R - y[i], -L - x[i])
        theta3[i] = np.arctan2(-R - y[i], -L - x[i])
        theta4[i] = np.arctan2( -R - y[i], L - x[i])
        
    return x, y, alpha

#x1, y1, alpha1 = collision_data_nospin(50, L, 0.5, 0.5)
x_spin, y_spin, alpha_spin = collision_data(50, L, 0.5, 0.5, 0, 1/2)

svals = np.linspace(-1, 1, 1000)
'''
plt.scatter(x1,y1)
plt.plot(x1,y1, c='green')
'''
plt.scatter(x_spin,y_spin, c='purple')
plt.plot(x_spin,y_spin, c='red')
plt.plot(L*svals, np.ones(len(svals)), c='black')
plt.plot(L*svals, -np.ones(len(svals)), c='black')
plt.plot(-np.cos(np.pi*svals/2) - L,np.sin(np.pi*svals/2), c='black')
plt.plot(np.cos(np.pi*svals/2) + L,-np.sin(np.pi*svals/2), c='black')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()

import numpy as np


def generate_collisions_spin(x0_vec, N, a = 1, r=1, MI_coeff=0.5):
    # x0_vec = [x0, y0, vx0, vy0, u0]
    x0, y0, vx0, vy0, u0 = x0_vec
    N = N+1
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
    alpha[0] = np.arctan2(vy0, vx0)
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(r - y[0], a - x[0])
    theta2[0] = np.arctan2(r - y[0], -a - x[0])
    theta3[0] = np.arctan2(-r - y[0], -a - x[0])
    theta4[0] = np.arctan2( -r - y[0], a - x[0])
    
    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (r - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = r
            
            vT = -vx[i - 1]
            vn = -vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = -vparr
            vy[i] = -vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] + a)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] + a)**2 + y[i - 1]**2 - r**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] + a])
            n = np.array([-(x[i] + a),-y[i]])
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
            t[i] = (-r - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = -r
            
            vT = vx[i - 1]
            vn = vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vparr
            vy[i] = vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] - a)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] - a)**2 + y[i - 1]**2 - r**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] - a])
            n = np.array([-(x[i] - a),-y[i]])
            v = np.array([vx[i - 1], vy[i - 1]])
            
            vT = np.dot(v, T)
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -np.dot(v,n)
            
            vnew = vparr*T + vperp*n
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vnew[0]
            vy[i] = vnew[1]
            alpha[i] = np.arctan2(vy[i], vx[i])

        theta1[i] = np.arctan2(r - y[i], a - x[i])
        theta2[i] = np.arctan2(r - y[i], -a - x[i])
        theta3[i] = np.arctan2(-r - y[i], -a - x[i])
        theta4[i] = np.arctan2( -r - y[i], a - x[i])
    
    return np.transpose(np.array([x, y, vx, vy]))
    

def bounce(x0_vec, a=1, r=1, MI_coeff=0.5):
    x0, y0, vx0, vy0, u0 = x0_vec
    N = 2
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
    alpha[0] = np.arctan2(vy0, vx0)
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(r - y[0], a - x[0])
    theta2[0] = np.arctan2(r - y[0], -a - x[0])
    theta3[0] = np.arctan2(-r - y[0], -a - x[0])
    theta4[0] = np.arctan2( -r - y[0], a - x[0])
    
    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1])%(2*np.pi) < (theta2[i - 1] - theta1[i - 1])%(2*np.pi):
            t[i] = (r - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = r
            
            vT = -vx[i - 1]
            vn = -vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = -vparr
            vy[i] = -vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta2[i - 1])%(2*np.pi) < (theta3[i - 1] - theta2[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] + a)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] + a)**2 + y[i - 1]**2 - r**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] + a])
            n = np.array([-(x[i] + a),-y[i]])
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
            t[i] = (-r - y[i - 1])/vy[i - 1]
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = -r
            
            vT = vx[i - 1]
            vn = vy[i - 1]
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -vn
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vparr
            vy[i] = vperp
            alpha[i] = np.arctan2(vy[i],vx[i])
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            b[i] = ((x[i - 1] - a)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] - a)**2 + y[i - 1]**2 - r**2)/(vx[i - 1]**2 + vy[i - 1]**2)
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*vx[i - 1]
            y[i] = y[i - 1] + t[i]*vy[i - 1]
            
            T = np.array([-y[i],x[i] - a])
            n = np.array([-(x[i] - a),-y[i]])
            v = np.array([vx[i - 1], vy[i - 1]])
            
            vT = np.dot(v, T)
            vparr = ((1 - MI_coeff)/(1 + MI_coeff))*vT - ((2*MI_coeff)/(1 + MI_coeff))*u[i - 1]
            vperp = -np.dot(v,n)
            
            vnew = vparr*T + vperp*n
            
            u[i] = -((1 - MI_coeff)/(1 + MI_coeff))*u[i - 1] - (2/(1 + MI_coeff))*vT
            vx[i] = vnew[0]
            vy[i] = vnew[1]
            alpha[i] = np.arctan2(vy[i], vx[i])

        theta1[i] = np.arctan2(r - y[i], a - x[i])
        theta2[i] = np.arctan2(r - y[i], -a - x[i])
        theta3[i] = np.arctan2(-r - y[i], -a - x[i])
        theta4[i] = np.arctan2( -r - y[i], a - x[i])
        
    return np.array([x[1], y[1], vx[1], vy[1], u[1]]), t[1]
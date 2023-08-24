import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# Width of rectangular part of the stadium:
L = 1

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
            b[i] = (x[i - 1] + L)*np.cos(alpha[i - 1]) + y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] + L)**2 + y[i - 1]**2 - 1
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] + L)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
            t[i] = (-1 - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = -1
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
            b[i] = (x[i - 1] - L)*np.cos(alpha[i - 1]) + y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] - L)**2 + y[i - 1]**2 - 1
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] - L)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
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
            b[i] = ((x[i - 1] + L)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] + L)**2 + y[i - 1]**2 - 1)/(vx[i - 1]**2 + vy[i - 1]**2)
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
            b[i] = ((x[i - 1] - L)*vx[i - 1] + y[i - 1]*vy[i - 1])/(vx[i - 1]**2 + vy[i - 1]**2)
            c[i] = ((x[i - 1] - L)**2 + y[i - 1]**2 - 1)/(vx[i - 1]**2 + vy[i - 1]**2)
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
        theta1[i] = np.arctan2(1 - y[i], L - x[i])
        theta2[i] = np.arctan2(1 - y[i], -L - x[i])
        theta3[i] = np.arctan2(-1 - y[i], -L - x[i])
        theta4[i] = np.arctan2( -1 - y[i], L - x[i])
        
    return x, y, u, vx, vy

# Define initial parameters
init_MI_coeff = 1/2
init_x = L
init_y = 0.5
init_theta = 0.5
init_u = 0
init_N = 50


x_spin, y_spin, u_spin, vx_spin, vy_spin = collision_data(init_N, init_x, init_y, init_theta, init_u, init_MI_coeff)

fig, ax = plt.subplots()
line, = ax.plot(x_spin, y_spin, lw=2, c='red')

svals = np.linspace(-1, 1, 1000)

plt.plot(L*svals, np.ones(len(svals)), c='black')
plt.plot(L*svals, -np.ones(len(svals)), c='black')
plt.plot(-np.cos(np.pi*svals/2) - L,np.sin(np.pi*svals/2), c='black')
plt.plot(np.cos(np.pi*svals/2) + L,-np.sin(np.pi*svals/2), c='black')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

axMI_coeff = fig.add_axes([0.25, 0.05, 0.65, 0.03])
MI_coeff_slider = Slider(
    ax=axMI_coeff,
    label='alpha',
    valmin=0,
    valmax=1,
    valinit=init_MI_coeff,
)

ax_x = fig.add_axes([0.25, 0.15, 0.65, 0.03])
x_slider = Slider(
    ax=ax_x,
    label='x0',
    valmin= -L - 1,
    valmax=L + 1,
    valinit=init_x,
)

ax_y = fig.add_axes([0.25, 0.1, 0.65, 0.03])
y_slider = Slider(
    ax=ax_y,
    label='y0',
    valmin=-1,
    valmax=1,
    valinit=init_y,
)

ax_theta = fig.add_axes([0.25, 0.2, 0.65, 0.03])
theta_slider = Slider(
    ax=ax_theta,
    label='theta0',
    valmin=0,
    valmax=2*np.pi,
    valinit=init_theta,
)

ax_u = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
u_slider = Slider(
    ax=ax_u,
    label="u",
    valmin=0,
    valmax=10,
    valinit=init_u,
    orientation="vertical"
)

ax_N = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
N_slider = Slider(
    ax=ax_N,
    label="N",
    valmin=1,
    valmax=100,
    valinit=init_N,
    orientation="vertical",
    valfmt='%0.0f'
)


# The function to be called anytime a slider's value changes
def update(val):
    x_spin, y_spin, u_spin, vx_spin, vy_spin = collision_data(int(N_slider.val), x_slider.val, y_slider.val, theta_slider.val, u_slider.val, MI_coeff_slider.val)
    line.set_xdata(x_spin)
    line.set_ydata(y_spin)
    fig.canvas.draw_idle()
    
MI_coeff_slider.on_changed(update)
x_slider.on_changed(update)
y_slider.on_changed(update)
theta_slider.on_changed(update)
u_slider.on_changed(update)
N_slider.on_changed(update)

resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    MI_coeff_slider.reset()
    x_slider.reset()
    y_slider.reset()
    theta_slider.reset()
    u_slider.reset()
    N_slider.reset()
button.on_clicked(reset)

plt.show()
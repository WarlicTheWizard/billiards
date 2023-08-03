import numpy as np
import matplotlib.pyplot as plt


def generate_collisions(x0, y0, alpha0, N_collisions, a=1, r=1, fname=None):
    # a is HALF the length of rectangular part of the stadium:
    # r is radius of semicircles
    # Number of collision points
    N = N_collisions + 1  # so that N_collisions is actual number of collisions
    t = np.zeros(N)  # times at collision
    x = np.zeros(N)  # x-values at collision
    y = np.zeros(N)  # y-values at collision
    phi = np.zeros(N)  # phi values at collision
    alpha = np.zeros(N)  # alpha values at collision
    b = np.zeros(N)  # b values at collision
    c = np.zeros(N)  # c values at collision
    theta1 = np.zeros(N)  # theta1 values at collision
    theta2 = np.zeros(N)  # theta2 values at collision
    theta3 = np.zeros(N)  # theta3 values at collision
    theta4 = np.zeros(N)  # theta4 values at collision

    # Initial values:
    t[0] = 0
    x[0] = x0  # 0.5
    y[0] = y0  # -0.5
    phi[0] = 0
    alpha[0] = alpha0  # np.pi/4 + 0.33
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(1 - y[0], a - x[0])
    theta2[0] = np.arctan2(1 - y[0], -a - x[0])
    theta3[0] = np.arctan2(-1 - y[0], -a - x[0])
    theta4[0] = np.arctan2(-1 - y[0], a - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1]) % (2*np.pi) < (theta2[i - 1] - theta1[i - 1]) % (2*np.pi):
            t[i] = (r - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = r
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta2[i - 1]) % (2*np.pi) < (theta3[i - 1] - theta2[i - 1]) % (2*np.pi):
            b[i] = (x[i - 1] + a)*np.cos(alpha[i - 1]) + \
                y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] + a)**2 + y[i - 1]**2 - r**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] + a)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        if (alpha[i - 1] - theta3[i - 1]) % (2*np.pi) < (theta4[i - 1] - theta3[i - 1]) % (2*np.pi):
            t[i] = (-r - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = -r
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta4[i - 1]) % (2*np.pi) < (theta1[i - 1] - theta4[i - 1]) % (2*np.pi):
            b[i] = (x[i - 1] - a)*np.cos(alpha[i - 1]) + \
                y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] - a)**2 + y[i - 1]**2 - r**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] - a)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        theta1[i] = np.arctan2(r - y[i], a - x[i])
        theta2[i] = np.arctan2(r - y[i], -a - x[i])
        theta3[i] = np.arctan2(-r - y[i], -a - x[i])
        theta4[i] = np.arctan2(-r - y[i], a - x[i])

    if fname:
        np.savetxt(fname, np.transpose(np.array([x, y, alpha])))

    return np.transpose(np.array([x, y, alpha]))


def plot_collisions(x, y, a=1, r=1):
    # x, y are coords of collision points
    svals = np.linspace(-1, 1, 1000)
    plt.scatter(x, y)
    plt.plot(x, y, c='green')  # for lines connecting
    # now plot stadium part
    # upper and lower lines
    plt.plot(a*svals, r*np.ones(len(svals)), c='black')
    plt.plot(a*svals, -r*np.ones(len(svals)), c='black')
    # semicircles
    plt.plot(-r*np.cos(np.pi*svals/2) - a, r*np.sin(np.pi*svals/2), c='black')
    plt.plot(r*np.cos(np.pi*svals/2) + a, -r*np.sin(np.pi*svals/2), c='black')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()


def bounce(x0, y0, alpha0, a=1, r=1):
    # just a single collision
    # use this function for calculating LCN
    N = 2
    t = np.zeros(N)  # times at collision
    x = np.zeros(N)  # x-values at collision
    y = np.zeros(N)  # y-values at collision
    phi = np.zeros(N)  # phi values at collision
    alpha = np.zeros(N)  # alpha values at collision
    b = np.zeros(N)  # b values at collision
    c = np.zeros(N)  # c values at collision
    theta1 = np.zeros(N)  # theta1 values at collision
    theta2 = np.zeros(N)  # theta2 values at collision
    theta3 = np.zeros(N)  # theta3 values at collision
    theta4 = np.zeros(N)  # theta4 values at collision

    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    phi[0] = 0
    alpha[0] = alpha0
    b[0] = 0
    c[0] = 0
    theta1[0] = np.arctan2(1 - y[0], a - x[0])
    theta2[0] = np.arctan2(1 - y[0], -a - x[0])
    theta3[0] = np.arctan2(-1 - y[0], -a - x[0])
    theta4[0] = np.arctan2(-1 - y[0], a - x[0])

    # Update formula:
    for i in range(1, N):
        if (alpha[i - 1] - theta1[i - 1]) % (2*np.pi) < (theta2[i - 1] - theta1[i - 1]) % (2*np.pi):
            t[i] = (r - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = r
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta2[i - 1]) % (2*np.pi) < (theta3[i - 1] - theta2[i - 1]) % (2*np.pi):
            b[i] = (x[i - 1] + a)*np.cos(alpha[i - 1]) + \
                y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] + a)**2 + y[i - 1]**2 - r**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] + a)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        if (alpha[i - 1] - theta3[i - 1]) % (2*np.pi) < (theta4[i - 1] - theta3[i - 1]) % (2*np.pi):
            t[i] = (-r - y[i - 1])/np.sin(alpha[i - 1])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = -r
            alpha[i] = -alpha[i - 1]
        if (alpha[i - 1] - theta4[i - 1]) % (2*np.pi) < (theta1[i - 1] - theta4[i - 1]) % (2*np.pi):
            b[i] = (x[i - 1] - a)*np.cos(alpha[i - 1]) + \
                y[i - 1]*np.sin(alpha[i - 1])
            c[i] = (x[i - 1] - a)**2 + y[i - 1]**2 - r**2
            t[i] = -b[i] + np.sqrt(b[i]**2 - c[i])
            x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
            y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
            phi[i] = np.arctan2(y[i], x[i] - a)
            alpha[i] = 2*phi[i] - alpha[i - 1] + np.pi
        theta1[i] = np.arctan2(r - y[i], a - x[i])
        theta2[i] = np.arctan2(r - y[i], -a - x[i])
        theta3[i] = np.arctan2(-r - y[i], -a - x[i])
        theta4[i] = np.arctan2(-r - y[i], a - x[i])

    return np.array([x[1], y[1], alpha[1]]), t[1]

import numpy as np
import matplotlib.pyplot as plt


def generate_collisions(x0, y0, alpha0, fname, N_collisions):

    # Number of collision points
    N = N_collisions

    t = np.zeros(N)  # times at collision
    x = np.zeros(N)  # x-values at collision
    y = np.zeros(N)  # y-values at collision
    phi = np.zeros(N)  # phi values at collision
    alpha = np.zeros(N)  # alpha values at collision
    b = np.zeros(N)  # b values at collision

    # Initial values:
    t[0] = 0
    x[0] = x0
    y[0] = y0
    phi[0] = 0
    alpha[0] = alpha0
    b[0] = x[0]*np.cos(alpha[0]) + y[0]*np.sin(alpha[0])
    c = x[0]**2 + y[0]**2 - 1  # c is only needed for one time-step

    # Update formulae for first collision
    t[1] = -b[0] + np.sqrt(b[0]**2 - c)
    x[1] = x[0] + t[1]*np.cos(alpha[0])
    y[1] = y[0] + t[1]*np.sin(alpha[0])
    phi[1] = np.arctan2(x[1], y[1])
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

    # May be uncommented to save collision data:
    np.savetxt(fname, np.transpose(np.array([x, y, alpha])))


# Plotting:
def plot_bounces(x, y):
    thetavals = np.linspace(0, 2*np.pi, 100)

    plt.scatter(x, y)
    plt.plot(np.cos(thetavals), np.sin(thetavals), color='black')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()

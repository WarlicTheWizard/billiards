import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def stadium_plot(x, y, a=1, r=1):
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


def rectangle_plot(x, y, a=1):
    plt.scatter(x, y)
    plt.plot(x, y, c='green')

    plt.gca().add_patch(Rectangle((-a, -1), 2*a, 2,
                        edgecolor='black',
                        facecolor='none'))
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()


def circle_plot(x, y, r=1):
    thetavals = np.linspace(0, 2*np.pi, 100)
    plt.scatter(x, y)
    plt.plot(x, y, c='green')
    plt.plot(np.cos(thetavals), np.sin(thetavals), color='black')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.show()

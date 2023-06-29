from stadium_generation import generate_collisions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from circle_animation_fixed import load_coordinates


dt = 0.005  # for the animation

def main():
    # explanation: cluster animation calculates collisions for a billiard starting at 
    # (x0, y0), as well as a ring of N_orbits many billiards delta away from the centre billiard
    # all billiards have same alpha0 (initial angle of velocity)
    # calculates and displays for N_collisions
    # Rectangle left bottom corner is at (-L, -1), right top is (L, 1)
    # default L=1, can specify to be different in cluster_animate
    # Semicircle part has radius 1. 

    N_collisions = 15
    N_orbits = 10
    delta = 10**(-5)
    L = 1
    # want special orbits
    # this orbit does the 'diamond' shape
    theta = np.arctan2(1, L+1)
    x0 = -L-1 + 10*delta*np.cos(theta)
    y0 = 10*delta*np.sin(theta)
    # cluster_animate(x0, y0, theta, N_orbits, N_collisions, delta=0.01, save_animation=True)
    x0 = -L-1 + 10*delta
    y0 = 0
    # cluster_animate(x0, y0, 0, N_orbits, N_collisions, delta=delta, save_animation=True)
    
    # random start
    # multiplying by (1-delta) ensures no orbits are outside boundary
    rng = np.random.default_rng()  # random number from [0, 1)
    x = (rng.random() * 2*L -L + rng.random() * 2 - 1) * (1-delta)
    if x>-L and x < L:
        # within 'rectangle' section
        y = (rng.random()*2 -1) * (1-delta)
    else:
        ext = np.abs(x) - L # distance into a semi circle
        maxy = np.sqrt(1 - ext**2)  # by pythagoras
        y = np.sign(x) * (rng.random() * (2*maxy) - maxy) * (1-delta)
    
    alpha = rng.random() * 2 * np.pi
    cluster_animate(x, y, alpha, N_orbits, N_collisions, delta=delta, save_animation=True)

def cluster_animate(x0, y0, alpha0, N_orbits, N_collisions, L=1, delta=10**(-5), save_animation=False):
    # make circle of points around (x0, y0)
    # be sure that the starting point (x0, y0) is more than delta from the edges
    starting_xs, starting_ys = [x0], [y0]
    if N_orbits != 0:
        angle_increment = 2 * np.pi / N_orbits
    for i in range(N_orbits):
        starting_xs.append(x0 + delta * np.cos(angle_increment * i))
        starting_ys.append(y0 + delta * np.sin(angle_increment * i))
    
    # now calculate collisions for each starting points, all with same angle
    cluster_coordinates, cluster_frames = [], []
    for i in range(N_orbits + 1):
        collisions = generate_collisions(starting_xs[i], starting_ys[i], alpha0, N_collisions, L=L)
        coordinates, n_frames = load_coordinates(collisions)
        cluster_coordinates.append(coordinates)
        cluster_frames.append(n_frames)
    
    
    fig = plt.figure()
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(xlim=(-L - 1.01, L + 1.01), ylim=(-1.01, 1.01))
    plt.axis('off')  # hide borders and axes
    ax.set_aspect(1)
    # draw background stadium
    svals = np.linspace(-1, 1, 1000)
    plt.plot(L*svals, np.ones(len(svals)), c='black')
    plt.plot(L*svals, -np.ones(len(svals)), c='black')
    plt.plot(-np.cos(np.pi*svals/2) - L,np.sin(np.pi*svals/2), c='black')
    plt.plot(np.cos(np.pi*svals/2) + L,-np.sin(np.pi*svals/2), c='black')
    ax.set_aspect('equal', adjustable='box')


    center, = plt.plot([], [], 'o', markersize=4, color='red')
    graph, = plt.plot([], [], 'o', markersize=2, color='blue')
    
    
    def init():
        graph.set_data([], [])
        center.set_data([], [])
        return graph, center
    
    def animate(i):
        # animate center particle
        x, y = cluster_coordinates[0][i]
        center.set_data([x], [y])
        x, y = [], []
        for coordinates in cluster_coordinates[1:]:
            xi, yi = coordinates[i]
            x.append(xi)
            y.append(yi)
        
        graph.set_data(x, y)
        return graph,
    
    n_frames = np.min(cluster_frames)  # some might have more frames than others
    ani = FuncAnimation(fig, animate, frames=n_frames, interval=1, blit=True)
    plt.show()
    if save_animation:
        ani.save('stadium_{}_{}_{}.mp4'.format(round(x0, 4), round(y0, 4), round(alpha0, 4)), fps=100)
        
if __name__=='__main__':
    main()
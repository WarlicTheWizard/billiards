from rectangle_generation import generate_collisions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

from circle_animation_fixed import load_coordinates


dt = 0.005

def main():
    # explanation: cluster animation calculates collisions for a billiard starting at 
    # (x0, y0), as well as a ring of N_orbits many billiards delta away from the centre billiard
    # all billiards have same alpha0 (initial angle of velocity)
    # calculates and displays for N_collisions
    # Length of rectangle is 2L, with bottom left corner at (-L, -1) 
    # and top right at (L, 1). Default L is 2
    
    N_collisions = 10
    N_orbits = 10
    delta = 0.01
    cluster_animate(0, 0, np.pi/4, N_orbits, N_collisions, delta=0.01, save_animation=True)
    cluster_animate(0.5, 0, np.pi/2, N_orbits, N_collisions, delta=delta, save_animation=True)
    rng = np.random.default_rng()  # random number from [0, 1)
    x = (rng.random() * 4 -2) * (1-delta)
    y = (rng.random() * 2 -1) * (1-delta)  # makes sure starting position is more than delta from boundary
    alpha = rng.random() * 2 * np.pi
    cluster_animate(x, y, alpha, N_orbits, N_collisions, delta=delta, save_animation=True)

def cluster_animate(x0, y0, alpha0, N_orbits, N_collisions, L=2, delta=10**(-5), save_animation=False):
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
    ax = plt.axes(xlim=(-L - 0.1, L + 0.1), ylim=(-1.1, 1.01))
    plt.axis('off')  # hide borders and axes
    ax.set_aspect(1)
    ax.add_patch(Rectangle((-L,-1),2*L,2,
                        edgecolor='black',
                        facecolor='none'))
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
        ani.save('rectangle_{}_{}_{}.mp4'.format(round(x0, 4), round(y0, 4), round(alpha0, 4)), fps=100)
        
if __name__=='__main__':
    main()
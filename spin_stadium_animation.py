from spin_stadium_generation import generate_collisions_spin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


dt = 0.005  # for the animation

def main():
    # explanation: cluster animation calculates collisions for a billiard starting at 
    # (x0, y0), as well as a ring of N_orbits many billiards delta away from the centre billiard
    # all billiards have same alpha0 (initial angle of velocity)
    # calculates and displays for N_collisions
    # Rectangle left bottom corner is at (-a, -r), right top is (a, r)
    # default a=1, can specify to be different in cluster_animate
    # Semicircle part has radius r=1. 

    N_collisions = 50
    N_orbits = 10
    delta = 10**(-5)
    a = 1
    # want special orbits
    # this orbit does the 'diamond' shape
    theta = np.arctan2(1, a+1)
    x0 = -a-1 + 10*delta*np.cos(theta)
    y0 = 10*delta*np.sin(theta)
    # cluster_animate(x0, y0, theta, N_orbits, N_collisions, delta=0.01, save_animation=True)
    x0 = -a-1 + 10*delta
    y0 = 0
    cluster_animate(x0, y0, 0, N_orbits, N_collisions, delta=delta, save_animation=True)
    
    # random start
    # multiplying by (1-delta) ensures no orbits are outside boundary
    rng = np.random.default_rng()  # random number from [0, 1)
    x = (rng.random() * 2*a -a + rng.random() * 2 - 1) * (1-delta)
    if x>-a and x < a:
        # within 'rectangle' section
        y = (rng.random()*2 -1) * (1-delta)
    else:
        ext = np.abs(x) - a # distance into a semi circle
        maxy = np.sqrt(1 - ext**2)  # by pythagoras
        y = np.sign(x) * (rng.random() * (2*maxy) - maxy) * (1-delta)
    
    alpha = rng.random() * 2 * np.pi
    # cluster_animate(x, y, alpha, N_orbits, N_collisions, delta=delta, save_animation=True)


def cluster_animate(x0, y0, alpha0, N_orbits, N_collisions, a=1, r=1, delta=10**(-5), save_animation=False, MI_coeff=0.5):
    # make circle of points around (x0, y0)
    # be sure that the starting point (x0, y0) is more than delta from the edges
    starting_xs, starting_ys = [x0], [y0]
    if N_orbits != 0:
        angle_increment = 2 * np.pi / N_orbits
    for i in range(N_orbits):
        starting_xs.append(x0 + delta * np.cos(angle_increment * i))
        starting_ys.append(y0 + delta * np.sin(angle_increment * i))
    
    # now calculate collisions for each starting points, all with same velocity
    # and same starting spin
    cluster_coordinates, cluster_frames = [], []
    for i in range(N_orbits + 1):
        x0_vec = [starting_xs[i], starting_ys[i], np.cos(alpha0), np.sin(alpha0), 1]
        collisions = generate_collisions_spin(x0_vec, N_collisions, a=a, r=r, MI_coeff=MI_coeff)
        coordinates, n_frames = load_coordinates(collisions)
        cluster_coordinates.append(coordinates)
        cluster_frames.append(n_frames)
    
    
    fig = plt.figure()
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(xlim=(-a - 1.01, a + 1.01), ylim=(-1.01, 1.01))
    plt.axis('off')  # hide borders and axes
    ax.set_aspect(1)
    # draw background stadium
    svals = np.linspace(-1, 1, 1000)
    plt.plot(a*svals, np.ones(len(svals)), c='black')
    plt.plot(a*svals, -np.ones(len(svals)), c='black')
    plt.plot(-np.cos(np.pi*svals/2) - a,np.sin(np.pi*svals/2), c='black')
    plt.plot(np.cos(np.pi*svals/2) + a,-np.sin(np.pi*svals/2), c='black')
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
        ani.save('stadium_{}_{}_{}.mp4'.format(round(x0, 4), round(y0, 4), round(alpha0, 4)), fps=200)
    
    
def load_coordinates(collisions_array):
    xedges = collisions_array[:, 0]
    yedges = collisions_array[:, 1]
    vxedges = collisions_array[:, 2]
    vyedges = collisions_array[:, 3]
    n_frames = 0
    x, y = xedges[0], yedges[0]
    initial_pos = np.array([x, y])
    coordinates = []
    coordinates.append(initial_pos)
    for i in range(len(xedges) - 1):
        x, y = xedges[i], yedges[i]
        initial_pos = np.array([x, y])
        v = np.array([vxedges[i], vyedges[i]])
        xnext, ynext = xedges[i+1], yedges[i+1]
        final_pos = np.array([xnext, ynext])
        max_length = np.linalg.norm(initial_pos - final_pos)
        v_mag = np.linalg.norm(v)
        d = 0
        while d < max_length:
            coordinates.append(coordinates[-1] + v*dt)
            d += v_mag * dt
            n_frames += 1
    return coordinates, n_frames



if __name__=='__main__':
    main()
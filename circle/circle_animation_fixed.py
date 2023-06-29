from collision_generation import generate_collisions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


dt = 0.005


def load_coordinates(fname):
    array = np.loadtxt(fname)
    xedges = array[:, 0]
    yedges = array[:, 1]
    alphaedges = array[:, 2]
    n_frames = 0
    x, y = xedges[0], yedges[0]
    initial_pos = np.array([x, y])
    coordinates = []
    coordinates.append(initial_pos)
    for i in range(len(xedges) - 1):
        x, y = xedges[i], yedges[i]
        initial_pos = np.array([x, y])
        alpha = alphaedges[i]
        v = np.array([np.cos(alpha), np.sin(alpha)])
        xnext, ynext = xedges[i+1], yedges[i+1]
        final_pos = np.array([xnext, ynext])
        max_length = np.linalg.norm(initial_pos - final_pos)
        t = 0
        while t < max_length:
            coordinates.append(coordinates[-1] + v*dt)
            t += dt
            n_frames += 1
    return coordinates, n_frames


x0, y0, alpha0 = 0.5, 0.5, np.pi/4
N_collisions = 10
epsilon = 10**-5  # initial separating distance
generate_collisions(x0, y0, alpha0, 'ball1.txt', N_collisions)
generate_collisions(x0+epsilon, y0, alpha0, 'ball2.txt', N_collisions)

ball1, frames1 = load_coordinates('ball1.txt')
ball2, frames2 = load_coordinates('ball2.txt')
n_frames = min(frames1, frames2)
distances = []
t = []
for i in range(n_frames):
    distances.append(np.linalg.norm(ball1[i] - ball2[i]))
    t.append(i * dt)

plt.plot(t, np.log(distances))
plt.show()

fig = plt.figure()
fig = plt.figure(figsize=(6, 6))
ax = plt.axes(xlim=(-1, 1), ylim=(-1, 1))
circ = plt.Circle((0, 0), 1, fill=False)
ax.set_aspect(1)
ax.add_artist(circ)
graph, = plt.plot([], [], 'o')


def animate(i):
    coord1, coord2 = ball1[i], ball2[i]
    x1, x2, y1, y2 = coord1[0], coord2[0], coord1[1], coord2[1]
    graph.set_data([x1, x2], [y1, y2])
    return graph


ani = FuncAnimation(fig, animate, frames=n_frames, interval=20)
#ani.save('bouncing.mp4', fps=50)
plt.show()

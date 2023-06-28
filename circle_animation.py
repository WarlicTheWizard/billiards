import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


array = np.loadtxt('edges.txt').T
xedges = array[:, 0]
yedges = array[:, 1]
alphaedges = array[:, 2]

# parameters
dt = 0.01
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


fig = plt.figure()
fig = plt.figure(figsize=(6, 6))
ax = plt.axes(xlim=(-1, 1), ylim=(-1, 1))
circ = plt.Circle((0, 0), 1, fill=False)
ax.set_aspect(1)
ax.add_artist(circ)
graph, = plt.plot([], [], 'o')


def animate(i):
    graph.set_data(coordinates[i])
    return graph


ani = FuncAnimation(fig, animate, frames=n_frames, interval=20)
plt.show()

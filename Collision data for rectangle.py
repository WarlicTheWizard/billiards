import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Width of rectangle:
L = 2

# Number of collision points
N = 50

t = np.zeros(N) # times at collision
x = np.zeros(N) # x-values at collision
y = np.zeros(N) # y-values at collision
alpha = np.zeros(N) # alpha values at collision
theta1 = np.zeros(N) # theta1 values at collision
theta2 = np.zeros(N) # theta2 values at collision
theta3 = np.zeros(N) # theta3 values at collision
theta4 = np.zeros(N) # theta4 values at collision

# Initial values:
t[0] = 0
x[0] = 0.5
y[0] = -0.5
alpha[0] = np.pi/4 + 0.33
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
        t[i] = (-L - x[i - 1])/np.cos(alpha[i - 1])
        x[i] = -L
        y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
        alpha[i] = np.pi - alpha[i - 1]
    if (alpha[i - 1] - theta3[i - 1])%(2*np.pi) < (theta4[i - 1] - theta3[i - 1])%(2*np.pi):
        t[i] = (-1 - y[i - 1])/np.sin(alpha[i - 1])
        x[i] = x[i - 1] + t[i]*np.cos(alpha[i - 1])
        y[i] = -1
        alpha[i] = -alpha[i - 1]
    if (alpha[i - 1] - theta4[i - 1])%(2*np.pi) < (theta1[i - 1] - theta4[i - 1])%(2*np.pi):
        t[i] = (L - x[i - 1])/np.cos(alpha[i - 1])
        x[i] = L
        y[i] = y[i - 1] + t[i]*np.sin(alpha[i - 1])
        alpha[i] = np.pi - alpha[i - 1]
    theta1[i] = np.arctan2(1 - y[i], L - x[i])
    theta2[i] = np.arctan2(1 - y[i], -L - x[i])
    theta3[i] = np.arctan2(-1 - y[i], -L - x[i])
    theta4[i] = np.arctan2( -1 - y[i], L - x[i])
    
# May be uncommented to save collision data:
#np.savetxt('rectangle_edges.txt',np.transpose(np.array([x,y,alpha])))

plt.scatter(x,y)
plt.plot(x,y, c='green')
plt.gca().add_patch(Rectangle((-L,-1),2*L,2,
                    edgecolor='black',
                    facecolor='none'))
plt.show()

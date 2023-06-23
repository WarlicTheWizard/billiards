import numpy as np
import matplotlib.pyplot as plt

# Width of rectangular part of the stadium:
L = 1

# Number of collision points
N = 50

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
x[0] = 0.5
y[0] = -0.5
phi[0] = 0
alpha[0] = np.pi/4 + 0.33
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

svals = np.linspace(-1, 1, 1000)

plt.scatter(x,y)
plt.plot(x,y, c='green')
plt.plot(L*svals, np.ones(len(svals)), c='black')
plt.plot(L*svals, -np.ones(len(svals)), c='black')
plt.plot(-np.cos(np.pi*svals/2) - L,np.sin(np.pi*svals/2), c='black')
plt.plot(np.cos(np.pi*svals/2) + L,-np.sin(np.pi*svals/2), c='black')
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()
import numpy as np
from stadium_generation import bounce
import matplotlib.pyplot as plt

e_mag = 1.9 * 10**(-7)  # length to which we normalize


def main():
    # starting x and theta
    # k, t = random_starts_LCN(100, 1000, a=1)
    LCN_vs_gamma(100, 100, 100)


def LCN_vs_gamma(N_gammas, N_starts, N_iters):
    # plot LCN as a function of gamma = a/r
    # N_starts: number of different initial positions to calculate LCN from
    # N_iters - how many collisions to go through for calculating LCN
    G = 1  # normalized area of stadium
    a_values = np.linspace(0, 1.7, N_gammas)
    # r as function of a and G
    r_values = (-4*a_values + np.sqrt(16*a_values**2+4*np.pi*G))/(2*np.pi)
    gamma_values = a_values/r_values
    LCNs = np.copy(a_values)
    for i in range(len(LCNs)):
        LCNs[i] = random_starts_LCN(
            N_starts, N_iters, a=a_values[i], r=r_values[i])
    plt.plot(gamma_values, LCNs)
    plt.title(r'Specular stadium billiard - LCN as function of $\gamma$ ratio')
    plt.xlabel(r'$\gamma$ ($ = a/r$, with area normalized)')
    plt.ylabel(r'LCN $\lambda$')


def random_starts_LCN(N_starts, N_iters, a=1, r=1):
    # calculate LCN many times for random starts within the stadium
    # return mean of LCNs calculated
    # random start
    # multiplying by (1-delta) ensures no orbits are outside boundary
    delta = 10**(-4)
    rng = np.random.default_rng()  # random number from [0, 1)
    x_array = []  # all starting vectors
    # generate random position in stadium

    for n in range(N_starts):
        x = (rng.random() * 2*a - a + rng.random() * 2*r - r) * (1-delta)
        if x > -a and x < a:
            # within 'rectangle' section
            y = (rng.random()*2*r - r) * (1-delta)
        else:
            ext = np.abs(x) - a  # distance into a semi circle
            maxy = np.sqrt(np.abs(r**2 - ext**2))  # by pythagoras
            y = np.sign(x) * (rng.random() * (2*maxy) - maxy) * (1-delta)
        alpha = rng.random() * 2 * np.pi
        x = np.array([x, y, alpha])
        x_array.append(x)

    LCN_array = []

    for x in x_array:
        LCN = calc_LCN(x, N_iters, a=a)
        LCN_array.append(LCN)

    return np.mean(LCN_array)


def calc_LCN(x, N_iters, a=1, r=1):
    # x is [x0, y0, theta0], as numpy array
    # generate random dx
    dx = np.random.rand(2) - 0.5
    dtheta = np.random.rand() * 2 * np.pi
    # normalize to make perturbation to specified small size
    dx = np.array([dx[0], dx[1], dtheta])
    dx = dx * e_mag / np.linalg.norm(dx)
    y = x + dx

    beta_array = np.ones(N_iters)
    t_array = []

    # we have analytical result for mean time
    mean_time = np.pi * (4*a*r + np.pi*r**2)/(2*np.pi*r + 4*a)

    for n in range(N_iters):

        x, y, t_x, t_y = bounce_endpoints(x, y, a=a)

        beta = e_mag / np.linalg.norm(x-y)  # calculate reduction factor
        # do reduction replacement
        dx = y - x
        y = x + dx * beta
        beta_array[n] = beta
        t_array.append([t_x, t_y])

    # calculate lyapunov exponent
    LCN = - np.sum(np.log(beta_array)) / (mean_time * N_iters)

    #t_array = np.array(t_array)
    #actual_mean = np.mean((t_array))
    #mean_diff = np.mean(np.abs(t_array[:,0] - t_array[:,1]))

    #print('actual mean time is ', actual_mean)
    #print('mean diff is', mean_diff)
    #print('analytical mean time is', mean_time)

    return LCN


def bounce_endpoints(x, y, a=1, r=1):
    # evolve the endpoints one collisions each
    x, t_x = bounce(*x, a=a, r=r)
    y, t_y = bounce(*y, a=a, r=r)

    return x, y, t_x, t_y


if __name__ == '__main__':
    main()

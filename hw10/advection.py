import numpy as np
import matplotlib.pyplot as plt

def myf(x):
    fac = 20
    return np.exp(-fac * (np.mod(x, 25) - 5)**2)

def advection(method=8):
    # Set up parameters
    a = np.sqrt(2)
    h = 0.05
    x = np.arange(0, 25, h)
    n = len(x)
    lambda_ = 0.8
    k = lambda_ * h / a
    tmax = 25
    u = myf(x)
    uold = np.copy(u)
    t = 0

    # Setup figure
    plt.ion()
    fig, ax = plt.subplots(figsize=(10, 5))
    exact_line, = ax.plot(x, u, 'r-', linewidth=2, label='Exact')
    num_line, = ax.plot(x, u, 'k-', linewidth=2, label='Numerical')
    ax.set_xlim([0, 25])
    ax.set_ylim([-0.5, 1.5])
    ax.set_title(method_name(method), fontsize=20)
    ax.legend()
    ax.grid(True)

    while t < tmax:
        ujp1 = np.roll(u, -1)
        ujm1 = np.roll(u, 1)

        if method == 1:  # Central difference
            unew = u - 0.5 * lambda_ * (ujp1 - ujm1)

        elif method == 2:  # Lax-Friedrichs
            unew = 0.5 * (ujm1 + ujp1 - lambda_ * (ujp1 - ujm1))

        elif method == 3:  # Upwind left
            unew = u - lambda_ * (u - ujm1)

        elif method == 4:  # Upwind right
            unew = u - lambda_ * (ujp1 - u)

        elif method == 5:  # Lax-Wendroff
            unew = u - 0.5 * lambda_ * (ujp1 - ujm1) + 0.5 * lambda_**2 * (ujp1 - 2*u + ujm1)

        elif method == 6:  # Beam-Warming left
            ujm2 = np.roll(u, 2)
            unew = u - 0.5 * lambda_ * (3*u - 4*ujm1 + ujm2) + 0.5 * lambda_**2 * (u - 2*ujm1 + ujm2)

        elif method == 7:  # Beam-Warming right
            ujp2 = np.roll(u, -2)
            unew = u - 0.5 * lambda_ * (-3*u + 4*ujp1 - ujp2) + 0.5 * lambda_**2 * (u - 2*ujp1 + ujp2)

        elif method == 8:  # Leap-frog
            if t < k:
                unew = u - 0.5 * lambda_ * (ujp1 - ujm1 + lambda_ * (ujp1 - 2*u + ujm1))
            else:
                unew = uold - lambda_ * (ujp1 - ujm1)
            uold = np.copy(u)

        # Update time and solution
        t += k
        u = unew

        # Update plot
        exact_line.set_ydata(myf(x - a*t))
        num_line.set_ydata(u)
        plt.pause(0.01)

    plt.ioff()
    plt.show()

def method_name(method):
    return {
        1: 'Central difference',
        2: 'Lax-Friedrichs',
        3: 'Upwind, left',
        4: 'Upwind, right',
        5: 'Lax-Wendroff',
        6: 'Beam-Warming, left',
        7: 'Beam-Warming, right',
        8: 'Leap-frog'
    }.get(method, 'Unknown method')

if __name__ == "__main__":
    advection(method=4)  # Change method number here

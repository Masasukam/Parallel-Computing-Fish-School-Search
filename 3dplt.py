import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = np.linspace(0, 10, 100)
y = np.linspace(0, 10, 100)
X, Y = np.meshgrid(x, y)

# a function that changes with time t
def f(X, Y, t):
    # Z = np.sin(0.01 * X * X + 0.05 * Y * Y + np.sin(t / 100) + 0.1 * Y + 0.5 * X)
    Z = np.sin(0.01 * X * X + 0.005 * Y * Y - 0.05 * X + 2 * np.sin(t / 100))
    return Z

steps = np.linspace(1, 1000, 50)  # arg3 times in the range [arg1, arg2]

for t in steps:
    Z = f(X, Y, t)  # Calculate Z for current time
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('f(x, y)')
    ax.set_title('3D Plot at t = {:.2f}'.format(t))

    plt.show()
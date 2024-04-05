import numpy as np
import matplotlib.pyplot as plt

# Example data
x = np.linspace(0, 10, 100)  # X-coordinates
y = np.linspace(0, 10, 100)  # Y-coordinates
X, Y = np.meshgrid(x, y)     # Create a grid of X and Y values
Z = np.sin(np.sqrt(X**2 + Y**2))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)')
plt.title('3D Plot of the function')

plt.show()

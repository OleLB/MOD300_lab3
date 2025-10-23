from matplotlib import pyplot as plt
import numpy as np


steps_per_walker = 1000000
walker_count = 5

# Generate random steps for each walker using numpy for efficiency
step_range = (-1, 1)
steps = np.random.uniform(step_range[0], step_range[1], (walker_count, steps_per_walker, 3))
# steps is an array of 3 dimensions: (walker_count, iterations, 3)

# Cumulative sum along each walker's steps to get positions over time, adds each 'point' to previous
positions = np.cumsum(steps, axis=1)

# Plotting the random walk
# Source for plot code: https://stackoverflow.com/questions/11541123/how-can-i-make-a-3d-line-plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

colors = ['r', 'g', 'b', 'y', 'c', 'm']

for i, w in enumerate(positions):
    x_vals = w[:, 0]
    y_vals = w[:, 1]
    z_vals = w[:, 2]

    ax.plot(x_vals, y_vals, z_vals, color=colors[i % len(colors)], label=f'Walker {i+1}')

ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

plt.legend()
plt.show()
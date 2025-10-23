import random
from matplotlib import pyplot as plt

class point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def set_point(self, x, y, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def get_point(self):
        return (self.x, self.y, self.z)

    def add_point(self, other_point):
        self.x += other_point.x
        self.y += other_point.y
        self.z += other_point.z


class walker:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.position = point(x, y, z)
        self.velocity = point()
        self.past_positions = []

    def move(self, point_change: point):
        self.past_positions.append(self.position.get_point())
        self.position.add_point(point_change)

    def get_position(self):
        return self.position.get_point()
    

def random_point(min_val=1, max_val=1):
    x = random.uniform(min_val, max_val)
    y = random.uniform(min_val, max_val)
    z = random.uniform(min_val, max_val)
    return point(x, y, z)


if __name__ == "__main__":
    iterations = 1000000
    walker_count = 5
    walkers = []

    for i in range(walker_count):
        w = walker()
        for j in range(iterations):
            step = random_point(-1, 1)
            w.move(step)
        walkers.append(w)

    # Plotting the random walk
    # Source for plot code: https://stackoverflow.com/questions/11541123/how-can-i-make-a-3d-line-plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colors = ['r', 'g', 'b', 'y', 'c', 'm']

    for i, w in enumerate(walkers):
        x_vals = [p[0] for p in w.past_positions]
        y_vals = [p[1] for p in w.past_positions]
        z_vals = [p[2] for p in w.past_positions]

        ax.plot(x_vals, y_vals, z_vals, color=colors[i % len(colors)], label=f'Walker {i+1}')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

    plt.legend()
    plt.show()
"""random walk, task 1 in topic 2, project 3 MOD300"""
import random
from matplotlib import pyplot as plt

# This is our implementation of the walkers without using numpy

class Point:
    """
    Point in 3D space
    """
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def set_point(self, x, y, z=0.0):
        """
        Set the point coordinates
        """
        self.x = x
        self.y = y
        self.z = z

    def get_point(self):
        """
        Get the point coordinates
        """
        return (self.x, self.y, self.z)

    def add_point(self, other_point):
        """
        Add another point to this point
        """
        self.x += other_point.x
        self.y += other_point.y
        self.z += other_point.z


class Walker:
    """
    A walker in 3D space
    """
    def __init__(self, x=0, y=0, z=0):
        self.position = Point(x, y, z)
        self.velocity = Point()
        self.past_positions = []

    def move(self, point_change: Point):
        """
        Move the walker by a point
        """
        self.past_positions.append(self.position.get_point())
        self.position.add_point(point_change)

    def get_position(self):
        """
        Get the position of the walker
        """
        return self.position.get_point()


def random_point(min_val=1, max_val=1):
    """
    Generate a random point in 3D space within given range
    """
    x = random.uniform(min_val, max_val)
    y = random.uniform(min_val, max_val)
    z = random.uniform(min_val, max_val)
    return Point(x, y, z)


def gen_walkers(iterations, walker_count, start_min=-100, start_max=100):
    """
    Generate a list of walkers and move them randomly
    """
    walkers_list = []
    for indx in range(walker_count):
        start = random_point(start_min, start_max)
        walker = Walker(start.x, start.y, start.z)
        # testing random starting points
        print(f"Walker {indx+1} start point: {start.x}, {start.y}, {start.z}")
        for _ in range(iterations):
            step = random_point(-1, 1)
            walker.move(step)
        walkers_list.append(walker)
    return walkers_list


if __name__ == "__main__":
    ITERATIONS = 10000
    WALKER_COUNT = 5

    walkers = gen_walkers(ITERATIONS, WALKER_COUNT)

    # Plotting the random walk
    # Source for plot code:
    # https://stackoverflow.com/questions/11541123/how-can-i-make-a-3d-line-plot
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

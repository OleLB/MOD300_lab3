"""Walker strategy sphere test file"""
import random
import matplotlib.pyplot as plt
import numpy
from main import (gen_spheres, gen_simulation_box, point_in_spheres,
    gen_random_points, get_atoms, get_dna_box, Point)

class Walker:
    """Walker object that moves in 3D space."""
    def __init__(self, box, x=0.0, y=0.0, z=0.0):
        self.position = Point(x, y, z)
        self.past_positions = []
        self.box = box

    def move(self, point_change):
        """
        Move the walker by a point, with wrapping
        """
        self.past_positions.append(self.position)
        new_pos = Point(
            self.position.x + point_change.x,
            self.position.y + point_change.y,
            self.position.z + point_change.z
        )

        # Handle wrapping around the box boundaries
        if new_pos.x > self.box.max_x:
            new_pos.x = self.box.min_x + (new_pos.x % self.box.max_x)
        elif new_pos.x < self.box.min_x:
            new_pos.x = self.box.max_x - (self.box.min_x - new_pos.x) % self.box.max_x

        if new_pos.y > self.box.max_y:
            new_pos.y = self.box.min_y + (new_pos.y % self.box.max_y)
        elif new_pos.y < self.box.min_y:
            new_pos.y = self.box.max_y - (self.box.min_y - new_pos.y) % self.box.max_y

        if new_pos.z > self.box.max_z:
            new_pos.z = self.box.min_z + (new_pos.z % self.box.max_z)
        elif new_pos.z < self.box.min_z:
            new_pos.z = self.box.max_z - (self.box.min_z - new_pos.z) % self.box.max_z

        self.position = new_pos

    def get_position(self):
        """
        Get the position of the walker
        """
        return self.position


def random_point(min_val=1, max_val=1):
    """
    Generate a random point in 3D space
    """
    x_random = random.uniform(min_val, max_val)
    y_random = random.uniform(min_val, max_val)
    z_random = random.uniform(min_val, max_val)
    return Point(x_random, y_random, z_random)


def gen_walkers(iterations, walker_count, walker_box):
    """
    Generate a list of walkers and move them randomly
    """
    walkers_list = []
    for _ in range(walker_count):
        walker = Walker(walker_box)
        starting_point = gen_random_points(walker_box, 1)[0]
        walker.position = starting_point
        for _ in range(iterations):
            step = random_point(-1, 1)
            walker.move(step)
        walkers_list.append(walker)

    return walkers_list


def accessible_volume_simulation(walker_count, steps_per_walker, simulation_box, spheres):
    """Run the accessible volume simulation."""

    w = gen_walkers(steps_per_walker, walker_count, simulation_box)
    box_volume = simulation_box.volume()
    points_in_spheres = 0

    for walkers in w:
        for past_position in walkers.past_positions:
            points_in_spheres += point_in_spheres(spheres, past_position)

    total_points = walker_count * steps_per_walker
    fraction = points_in_spheres / total_points
    sphere_volume_estimate = fraction * box_volume
    print(f"Estimated volume of spheres: {sphere_volume_estimate}")

    actual_sphere_volume = 0
    for sphere in spheres:
        actual_sphere_volume += sphere.volume()

    print(f"Actual volume of spheres: {actual_sphere_volume}")
    print(f"Box volume: {box_volume}")
    accessible_volume_estimate = box_volume - sphere_volume_estimate
    print("Estimated accessible volume: "
          + f"(box volume - estimated volume of spheres) = {accessible_volume_estimate}")

    error = actual_sphere_volume - sphere_volume_estimate

    return sphere_volume_estimate, error

    # Source code for 3d subplots:
    # https://matplotlib.org/stable/gallery/mplot3d/subplot3d.html
    # fig1 = plt.figure()
    # ax = fig1.add_subplot(111, projection='3d')
    # ax.set_title('Simulation of the spheres (no walkers)')

    # # Taken from plot_points_in_spheres function in main.py
    # for sph in spheres:
    #     theta = numpy.linspace(0, 2 * numpy.pi, 100)
    #     phi = numpy.linspace(0, numpy.pi, 50)
    #     theta, phi = numpy.meshgrid(theta, phi)
    #     x_sph = sph.rad * numpy.sin(phi) * numpy.cos(theta) + sph.point.x
    #     y_sph = sph.rad * numpy.sin(phi) * numpy.sin(theta) + sph.point.y
    #     z_sph = sph.rad * numpy.cos(phi) + sph.point.z
    #     ax.plot_surface(x_sph, y_sph, z_sph, cmap='viridis', alpha=0.8)

    # ax.set_xlabel('X Axis')
    # ax.set_ylabel('Y Axis')
    # ax.set_zlabel('Z Axis')

    # colors = ['r', 'g', 'b', 'y', 'c', 'm']

    # fig2 = plt.figure()
    # ax = fig2.add_subplot(111, projection='3d')
    # ax.set_title('Walker paths')

    # for i, w in enumerate(w):
    #     x_vals = [p.x for p in w.past_positions]
    #     y_vals = [p.y for p in w.past_positions]
    #     z_vals = [p.z for p in w.past_positions]

    #     ax.plot(x_vals, y_vals, z_vals, color=colors[i % len(colors)])  # label=f'Walker {i+1}'

    # ax.set_xlabel('X Axis')
    # ax.set_ylabel('Y Axis')
    # ax.set_zlabel('Z Axis')

    # # plt.legend()
    # plt.show()


if __name__ == "__main__":

    WALKER_COUNT = 250
    STEPS = 100
    results = []

    # Spheres simulation
    # SIM_BOX = gen_simulation_box(50, 100, 50)
    # spheres = gen_spheres(SIM_BOX, 10)
    # for i in range(5):
    #     result, _ = accessible_volume_simulation(WALKER_COUNT, STEPS, SIM_BOX, spheres)
    #     results.append(result)
    # average = sum(results) / len(results)
    # print(f"Average result over 5 runs: {average}")

    # DNA simulation
    atoms = get_atoms()
    dna_box = get_dna_box(atoms)
    for i in range(5):
        result, _ = accessible_volume_simulation(WALKER_COUNT, STEPS, dna_box, atoms)
        results.append(result)
    average = sum(results) / len(results)
    print(f"Average result over 5 runs: {average}")

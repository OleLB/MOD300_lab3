"""Main .py file for MOD300 lab 3"""
# angstrom: a unit of length equal to one hundred-millionth of a centimetre, 10−10metre
# 1 angstrom is 100 pm (pictometer)
import random
import math
import numpy
from matplotlib import pyplot as plt

# atom_properties = {"atomic symbol: [atomic radius (pm), color in plot], ..."}
atom_properties = {
    "H":[120, "blue"], 
    "O":[152, "red"], 
    "P":[180, "yellow"], 
    "C":[170, "black"], 
    "N":[155, "green"]
}

class Point:
    """
    Point object in 3D
    """
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __sub__(self, other):
        """Ensure correct behavior when doing 'point - point' """
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

    def __str__(self):
        """Clean printing of point objects"""
        return f"x: {self.x}\ny: {self.y}\nz: {self.z}\n"

    # Source for equality and hashing:
    # https://stackoverflow.com/questions/1227121/compare-object-instances-for-equality-by-their-attributes
    def __eq__(self, other):
        """Check if two points are equal by comparing coordinates"""
        if not isinstance(other, Point):
            return False
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        """Allows Point objects to be used in sets and as dict keys"""
        return hash((self.x, self.y, self.z))


class Atom:
    """
    Atom object
    symbol: It's letter from the periodic table
    xyz: It's position
    """
    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.point = Point(x, y, z)
        self.rad = atom_properties[symbol][0]/100 # convert to angstrom units
        self.color = atom_properties[symbol][1]

    def volume(self) -> float:
        """Calculate the volume of the atom as a sphere"""
        volume = (4/3)*math.pi*self.rad**3
        return volume

    def __str__(self):
        return f'''
Atomic symbol: {self.symbol}
Atomic radius: {self.rad} angstrom
Position: 
{self.point}
'''


class Sphere:
    """
    Sphere object
    point: position of the sphere
    rad: radius of the sphere
    """
    def __init__(self, point, rad):
        self.point = point
        self.rad = rad

    def __str__(self):
        return f"Sphere radius: {self.rad}\nSphere location:\n{self.point}"

    def volume(self) -> float:
        """Calculate the volume of the sphere"""
        volume = (4/3)*math.pi*self.rad**3
        return volume


class SimulationBox:
    """
    Simulation box object
    Contains the coordinate limits for a box in 3D space
    """
    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z):
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z

    def volume(self):
        """Calculate the volume of the SimulationBox"""
        x_dist = self.max_x - self.min_x
        y_dist = self.max_y - self.min_y
        z_dist = self.max_z - self.min_z

        return x_dist*y_dist*z_dist

    def __str__(self):
        return f'''
min x: {self.min_x}
max x: {self.max_x}
min y: {self.min_y}
max y: {self.max_y}
min z: {self.min_z}
max z: {self.max_z}
'''


# Task 0
def gen_simulation_box(x, y, z):
    """
    Generates a dictionary that represents a simulation box with the sizes passed as arguments
    ---
    Returns a dictionary
    """
    return SimulationBox(0, x, 0, y, 0, z)



def gen_spheres(simulation_box, sphere_count) -> list[Sphere]:
    """
    Generates an ammount of spheres with random position and radius
    ---
    Returns a list of spheres
    """
    #! Currently it can generate overlapping spheres
    rand_spheres = []
    for _ in range(sphere_count):
        radius = random.uniform(1, 2)
        point = gen_random_points(simulation_box, 1, radius)[0]
        rand_spheres.append(Sphere(point, radius))

    return rand_spheres


# We learned how to calculate vector magnitude with numpy here:
# https://www.geeksforgeeks.org/python/how-to-get-the-magnitude-of-a-vector-in-numpy/
def point_in_spheres(spheres_or_atoms, point):
    """
    Calculate the distance from a point to every atom using numpy.
    If the magnitude of the vector between the point and any atom ...
      is smaller than the atomic radius of that atom, the point is in the DNA
    Vector between points = point1 - point2

    The first argument can be a list of sphere objects OR a list of atoms, both work
    ---
    Returns a bool
    """
    in_sphere_counter = 0
    # Turn to coordinates into numpy arrays
    dna_coords = numpy.array([[a.point.x, a.point.y, a.point.z] for a in spheres_or_atoms])
    single_point = numpy.array([point.x, point.y, point.z])

    vectors = dna_coords - single_point     # the single point gets subtracted from all dna points
    # print(f"vectors: {vectors[:5]}")

    magnitudes = numpy.linalg.norm(vectors, axis=1)

    # Compare each magnitude with the correct atomic radius
    for i,magnitude in enumerate(magnitudes):
        if magnitude <= spheres_or_atoms[i].rad:
            in_sphere_counter += 1
    return in_sphere_counter


def get_atoms() -> list:
    """
    Reads a file of DNA data and creates atom objects
    ---
    Returns a list of atom objects
    """
    atom_list = []

    with open("dna_coords.txt", encoding="utf-8") as coords:
        for line in coords:
            cur_atom_data = line.split()
            atom_list.append(Atom(*cur_atom_data))
    return atom_list


def get_dna_box(dna_atoms) -> SimulationBox:
    """
    Finds the lowest and highest xyz values from the atom coordinates.
    ---
    Returns a dictionary
    """
    smallest_x = math.inf
    largest_x = -math.inf
    smallest_y = math.inf
    largest_y = -math.inf
    smallest_z = math.inf
    largest_z = -math.inf

    for atom in dna_atoms:
        # Check x
        if atom.point.x > largest_x:
            largest_x = atom.point.x
        elif atom.point.x < smallest_x:
            smallest_x = atom.point.x

        # Check y
        if atom.point.y > largest_y:
            largest_y = atom.point.y
        elif atom.point.y < smallest_y:
            smallest_y = atom.point.y

        # Check z
        if atom.point.z > largest_z:
            largest_z = atom.point.z
        elif atom.point.z < smallest_z:
            smallest_z = atom.point.z

    # Make the box slightly larger to fit the atomic radius

    smallest_x = math.floor(smallest_x-5)
    largest_x = math.ceil(largest_x+5)
    smallest_y = math.floor(smallest_y-5)
    largest_y = math.ceil(largest_y+5)
    smallest_z = math.floor(smallest_z-5)
    largest_z = math.ceil(largest_z+5)

    return SimulationBox(smallest_x, largest_x, smallest_y, largest_y, smallest_z, largest_z)


def gen_random_points(box_points: SimulationBox, number_of_points: int, radius=0) -> list:
    """
    Generates points randomly placed within the simulation box.
    Can take the radius of a sphere to ensure that 
        no point of the sphere is outside the simulation box
    ---
    Return the points as a list.
    """
    points = []
    for _ in range(number_of_points):
        x = random.uniform(box_points.min_x+radius, box_points.max_x-radius)
        y = random.uniform(box_points.min_y+radius, box_points.max_y-radius)
        z = random.uniform(box_points.min_z+radius, box_points.max_z-radius)
        points.append(Point(x, y, z))

    return points


def plot_points_and_spheres(spheres_to_plot, points):
    """Function to make a plot of points and spheres"""
    # Source for sphere plotting code: https://likegeeks.com/3d-sphere-python/
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    theta = numpy.linspace(0, 2 * numpy.pi, 100)
    phi = numpy.linspace(0, numpy.pi, 50)
    theta, phi = numpy.meshgrid(theta, phi)
    for sph in spheres_to_plot:
        x = sph.rad * numpy.sin(phi) * numpy.cos(theta) + sph.point.x
        y = sph.rad * numpy.sin(phi) * numpy.sin(theta) + sph.point.y
        z = sph.rad * numpy.cos(phi) + sph.point.z
        ax.plot_surface(x, y, z, cmap='viridis', alpha=0.8)

    for point in points[:250]:
        ax.scatter(xs=point.x, ys=point.y, zs=point.z, color='blue')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

    # plt.legend()
    plt.show()



###########    Topic 2, task 1 and 2: Random Walkers     #################

### Task 1

class Walker:
    """Walker object that moves in 3D space."""
    def __init__(self, x=0.0, y=0.0, z=0.0, box=None, step_size=0):
        self.position = Point(x, y, z)
        self.past_positions = []
        self.box = box
        self.step_size = step_size

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

        if self.box:
            # Handle wrapping around the box boundaries if a box is defined
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

        # Adjusting to STEP_SIZE
        if self.step_size > 0:
            new_pos.x = round(new_pos.x / self.step_size) * self.step_size
            new_pos.y = round(new_pos.y / self.step_size) * self.step_size
            new_pos.z = round(new_pos.z / self.step_size) * self.step_size

            # rounding to avoid floating point precision issues
            new_pos = Point(
                round(new_pos.x, 1),
                round(new_pos.y, 1),
                round(new_pos.z, 1)
            )

        self.position = new_pos

    def set_new_position(self, point):
        """
        Set a new position for the walker
        """
        self.past_positions.append(self.position)
        self.position = point


def random_step(min_val=1, max_val=1):
    """
    Generate a random point in 3D space within given range
    """
    x = random.uniform(min_val, max_val)
    y = random.uniform(min_val, max_val)
    z = random.uniform(min_val, max_val)
    return Point(x, y, z)


def random_walk(iterations=10000, walker_count=5):
    """
    Generate a list of walkers and move them randomly
    The plot contains some long straight lines because the walker is sent to the 
        opposite side of the box when it hits a wall.
    """
    sim_box = gen_simulation_box(200, 200, 200)
    walkers_list = []
    for indx in range(walker_count):
        start = gen_random_points(sim_box, 1)[0]
        walker = Walker(start.x, start.y, start.z)
        # testing random starting points
        # print(f"Walker {indx+1} start point: {start.x}, {start.y}, {start.z}")
        for _ in range(iterations):
            step = random_step(-1, 1)
            walker.move(step)
        walkers_list.append(walker)
    
    # Plotting the random walk
    # Source for plot code:
    # https://stackoverflow.com/questions/11541123/how-can-i-make-a-3d-line-plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colors = ['r', 'g', 'b', 'y', 'c', 'm']

    for i, w in enumerate(walkers_list):
        x_vals = [p.x for p in w.past_positions]
        y_vals = [p.y for p in w.past_positions]
        z_vals = [p.z for p in w.past_positions]

        ax.plot(x_vals, y_vals, z_vals, color=colors[i % len(colors)], label=f'Walker {i+1}')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

    plt.legend()
    plt.show()


### Task 2

def random_walk_fast(steps_per_walker = 10000, walker_count = 5):
    """walker function with a focus on efficient use of numpy"""

    # Generate random steps for each walker using numpy for efficiency
    step_range = (-1, 1)
    steps = numpy.random.uniform(step_range[0], step_range[1], (walker_count, steps_per_walker, 3))
    # steps is an array of 3 dimensions: (walker_count, iterations, 3)

    # Cumulative sum along each walker's steps to get positions over time,
    #   adds each 'point' to previous
    positions = numpy.cumsum(steps, axis=1)

    # Plotting
    # Source for plot code:
    # https://stackoverflow.com/questions/11541123/how-can-i-make-a-3d-line-plot
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







###########    Topic 2, task 5: Surface area    #################

STEP_SIZE = 0.2  # Angstroms

def find_empty_position(box, spheres_to_avoid):
    """
    Find a random position in the box that is not inside any sphere.
    Ensure the position's x, y, z values are multiples of STEP_SIZE.
    """
    while True:
        point = gen_random_points(box, 1)[0]

        # Adjusting to STEP_SIZE
        point.x = round(point.x / STEP_SIZE) * STEP_SIZE
        point.y = round(point.y / STEP_SIZE) * STEP_SIZE
        point.z = round(point.z / STEP_SIZE) * STEP_SIZE

        # rounding to avoid floating point precision issues
        snapped_point = Point(
            round(point.x, 1),
            round(point.y, 1),
            round(point.z, 1)
        )
        if not point_in_spheres(spheres_to_avoid, snapped_point):
            # print("Found empty position:", snapped_point)
            return snapped_point


def spawn_walkers(count, box, spheres):
    """
    Spawn walkers at random empty positions in the box
    """
    walkers = []
    for _ in range(count):
        pos = find_empty_position(box, spheres)
        walker = Walker(pos.x, pos.y, pos.z, box, STEP_SIZE)
        walkers.append(walker)
    return walkers


def move_walker(walker, box, spheres, steps):
    """
    Move a walker randomly and check for surface hits
    """
    surface_blocks = set()
    move_options = [
        Point(STEP_SIZE, 0.0, 0.0), Point(-STEP_SIZE, 0.0, 0.0),
        Point(0.0, STEP_SIZE, 0.0), Point(0.0, -STEP_SIZE, 0.0),
        Point(0.0, 0.0, STEP_SIZE), Point(0.0, 0.0, -STEP_SIZE)
    ]
    for _ in range(steps):  
        move_choice = random.choice(move_options)
        walker.move(move_choice)

        # Check if the walker is inside any sphere
        if point_in_spheres(spheres, walker.position):
            surface_blocks.add(walker.position)
            # Move the walker to a new random position in the box
            walker.set_new_position(find_empty_position(box, spheres))
    return surface_blocks


def make_deterministic_spheres():
    """Create a set of deterministic spheres for testing."""
    spheres_list = []
    spheres_list.append(Sphere(Point(4.186925985859521, 4.15911461644369, 2.3552084360025116), 1.8887927910439066))
    spheres_list.append(Sphere(Point(5.683792634079552, 5.5375850544203065, 5.247502976020097), 1.5296427088887965))
    spheres_list.append(Sphere(Point(1.8917794544374336, 2.6501098749940146, 7.729438743292647), 1.813504604313081))
    return spheres_list


def estimate_surface_area(deterministic=True):
    """This function estimates the total surface area of spheres using walkers"""
    steps_per_walker = 3000
    walker_count = 300
    simulation_box = gen_simulation_box(10, 10, 10)

    if deterministic:
        spheres_list = make_deterministic_spheres()
    else:
        spheres_list = gen_spheres(simulation_box, 3)

    # Calculate actual surface area
    actual_surface_area = 0
    for sph in spheres_list:
        actual_surface_area += 4 * math.pi * (sph.rad ** 2)

    # Spawn walkers
    walker_list = spawn_walkers(walker_count, simulation_box, spheres_list)

    # Move walkers and record surface hits
    recorded_positions = set()
    for w in walker_list:
        surface_hits = move_walker(w, simulation_box, spheres_list, steps_per_walker)
        recorded_positions.update(surface_hits)

    # Get the past_positions of all walkers and remove duplicates
    all_visited_positions = set()
    for w in walker_list:
        all_visited_positions.update(w.past_positions)

    # Calculate estimated surface area
    volume_covered = len(all_visited_positions) * STEP_SIZE**3
    area_multiplier = 1 / STEP_SIZE
    fraction_of_box_covered = volume_covered / simulation_box.volume()
    fraction_of_surface_hits = len(recorded_positions) / len(all_visited_positions)

    points = simulation_box.volume() / (STEP_SIZE**3)
    estimated_surface_area_points = fraction_of_surface_hits * points
    estimated_surface_area = estimated_surface_area_points * (STEP_SIZE**2) * area_multiplier

    print("Fraction of box covered:", fraction_of_box_covered)
    print("Surface area hits recorded:", len(recorded_positions))
    print(f"Actual Surface Area: {actual_surface_area}")
    print(f"Estimated Surface Area: {estimated_surface_area}")

    # Plot the recorded positions
    xs = [p.x for p in recorded_positions]
    ys = [p.y for p in recorded_positions]
    zs = [p.z for p in recorded_positions]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs, c='b', marker='o', label='Surface Hits')

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_zlim(0, 10)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.legend()
    plt.show()





if __name__ == "__main__":

    estimate_surface_area()

    exit()

    # Verify that point subtraction works as expected
    print('-'*50+"Test point subtraction"+'-'*50)
    point1 = Point(5,2,1)
    point2 = Point(2,2,2)
    point3 = point1 - point2
    test_point = Point(5-2,2-2,1-2)
    assert (point3.x == test_point.x), "point subtraction failed"
    assert (point3.y == test_point.y), "point subtraction failed"
    assert (point3.z == test_point.z), "point subtraction failed"
    print('Subtraction assert passed')

    # Verify that simulation box is correctly created
    print('-'*50+"Test simulation box generation"+'-'*50)
    atoms = get_atoms()
    print(atoms[0])
    simulation_box_points = get_dna_box(atoms)
    #print(simulation_box_points)
    box_volume = simulation_box_points.volume()
    assert box_volume == 26970, "Simulation box is not correct"
    #print(f"Simulation box volume: {box_volume}")
    print("Simulation box assert passed")

    # Test point_in_spheres function
    print('-'*50+"Test point_in_spheres function"+'-'*50)
    fake_point = Point(1,1,1)  # not in dna
    real_point = Point(-48, 1.74, -1.22) # in dna
    assert not point_in_spheres(atoms, fake_point), "points_in_dna not working"
    assert point_in_spheres(atoms, real_point), "points_in_dna not working"
    print("point_in_spheres assert passed")


    # Test random point generation
    print('-'*50+"Test random point function"+'-'*50)
    NUMBER_OF_POINTS = 10
    random_points = gen_random_points(simulation_box_points, NUMBER_OF_POINTS)
    assert len(random_points) == NUMBER_OF_POINTS
    for p in random_points:
        assert p.x > simulation_box_points.min_x, "x coordinate not in simulation box"
        assert p.x < simulation_box_points.max_x, "x coordinate not in simulation box"
        assert p.y > simulation_box_points.min_y, "y coordinate not in simulation box"
        assert p.x < simulation_box_points.max_y, "y coordinate not in simulation box"
        assert p.z > simulation_box_points.min_z, "z coordinate not in simulation box"
        assert p.x < simulation_box_points.max_z, "z coordinate not in simulation box"
    print("Random coordinate assert passed")


    # Test gen spheres
    print('-'*50+"Test random sphere function"+'-'*50)
    spheres = gen_spheres(simulation_box_points, 3)
    for s in spheres:
        print(s)

    # plot_points_and_spheres(atoms[:20], [])

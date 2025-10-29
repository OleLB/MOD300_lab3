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
        """Calculate the volume of the sphere"""
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


if __name__ == "__main__":

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

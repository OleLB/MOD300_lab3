# angstrom: a unit of length equal to one hundred-millionth of a centimetre, 10−10metre
# 1 angstrom is 100 pm (pictometer)
import math

# atom_properties = {"atomic symbol: [atomic radius (pm), color in plot], ..."}
atom_properties = {"H":[120, "blue"], "O":[152, "red"], "P":[180, "yellow"], "C":[170, "black"], "N":[155, "green"]}

class point:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __sub__(self, other):
        """Ensure correct behavior when doing 'point - point' """
        return point(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __str__(self):
        """Clean printing of point objects"""
        return f"x: {self.x}\ny: {self.y}\nz: {self.z}\n"


class atom:
    def __init__(self, type, x, y, z):
        self.type = type
        self.coords = point(x, y, z)
        self.rad = atom_properties[type][0]/100 # convert to angstrom units
        self.color = atom_properties[type][1]


def get_atoms() -> list:
    atoms = []

    with open("dna_coords.txt") as coords:
        for line in coords:
            cur_atom_data = line.split()
            atoms.append(atom(*cur_atom_data))
    return atoms


def get_box_size(atoms) -> dict:
    smallest_x = math.inf
    largest_x = -math.inf
    smallest_y = math.inf
    largest_y = -math.inf
    smallest_z = math.inf
    largest_z = -math.inf

    for atom in atoms:
        # Check x
        if atom.coords.x > largest_x:
            largest_x = atom.coords.x
        elif atom.coords.x < smallest_x:
            smallest_x = atom.coords.x

        # Check y
        if atom.coords.y > largest_y:
            largest_y = atom.coords.y
        elif atom.coords.y < smallest_y:
            smallest_y = atom.coords.y

        # Check z
        if atom.coords.z > largest_z:
            largest_z = atom.coords.z
        elif atom.coords.z < smallest_z:
            smallest_z = atom.coords.z

    # Make the box slightly larger to fit the atomic radius
    box_points = {
        "smallest_x":math.floor(smallest_x-5),
        "largest_x":math.ceil(largest_x+5),
        "smallest_y":math.floor(smallest_y-5),
        "largest_y":math.ceil(largest_y+5),
        "smallest_z":math.floor(smallest_z-5),
        "largest_z":math.ceil(largest_z+5)
    }

    return box_points


def calc_box_volume(box_points) -> float:
    x_size = box_points["largest_x"] - box_points["smallest_x"]
    y_size = box_points["largest_y"] - box_points["smallest_y"]
    z_size = box_points["largest_z"] - box_points["smallest_z"]

    return x_size*y_size*z_size



if __name__ == "__main__":

    # Verify that point subtraction words as expected
    point1 = point(5,2,1)
    point2 = point(2,2,2)
    point3 = point1 - point2
    assert_point = point(5-2,2-2,1-2)
    
    assert (point3.x == assert_point.x and point3.y == assert_point.y and point3.z == assert_point.z), "point subtraction failed"

    atoms = get_atoms()
    box_points = get_box_size(atoms)
    print(box_points)
    print(calc_box_volume(box_points))
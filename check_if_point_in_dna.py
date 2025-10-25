# If the magnitude of the vector between the point and any 'atom' in the dna is smaller than the pm for that atom, the point is in the DNA
# Vector between points = point1 - point2

import numpy
from read_file import get_atoms, point


# We learned how to calculate vector magnitude with numpy here: https://www.geeksforgeeks.org/python/how-to-get-the-magnitude-of-a-vector-in-numpy/        
def point_in_dna(atoms, point):
    """Calculate the distance from a point to every atom using numpy"""
    # Turn to coordinates into numpy arrays
    dna_coords = numpy.array([[a.coords.x, a.coords.y, a.coords.z] for a in atoms])
    single_point = numpy.array([point.x, point.y, point.z])

    vectors = dna_coords - single_point     # the single point gets subtracted from all dna points
    # print(f"vectors: {vectors[:5]}")

    magnitudes = numpy.linalg.norm(vectors, axis=1)

    for i,magnitude in enumerate(magnitudes):    # This is to compare each magniture with the correct atomic radius
        if magnitude <= atoms[i].rad:
            return True
    return False


if __name__ == "__main__":
    atoms = get_atoms()
    fake_point = point(1,1,1)  # not in dna
    real_point = point(-48, 1.74, -1.22) # in dna
    assert False == point_in_dna(atoms, fake_point), "points_in_dna not working"
    assert True == point_in_dna(atoms, real_point), "points_in_dna not working"
""" Functions for testing """
import unittest

from main import Point, gen_random_points, get_dna_box, get_atoms, SimulationBox

class Tests(unittest.TestCase):
    """
        Class for testing using unittest.TestCase as an interface. 
    """

    def test_point_subtraction_1(self):
        """
            Test that point subtraction works.
        """
        p1 = Point(x=2, y=0, z=3)
        p2 = Point(x=4, y=2, z=3)
        result = p1 - p2

        self.assertEqual(result.x == -2, result.y == -2, result.z == 0)

    def test_assert_gen_random_point(self):
        """
            Test that assert number_of_points >= 0 works. 
        """
        sim_box = SimulationBox(0, 10, 0, 10, 0, 10)
        with self.assertRaises(AssertionError):
            gen_random_points(sim_box, -1, 2)

    def test_eq_gen_random_points(self):
        """
            Test that gen_random_points() returns the same amount 
            of points as parameter 'number_of_points'.
        """
        sim_box = SimulationBox(0, 10, 0, 10, 0, 10)
        self.assertEqual(10, len(gen_random_points(sim_box, 10)))

    def test_assert_get_dna_box(self):
        """
            Test that assert dna_atoms works. 
        """
        with self.assertRaises(AssertionError):
            get_dna_box(0)

    def test_get_atoms(self):
        """
            Test that get_atoms() returns a list.
        """
        self.assertIsInstance(get_atoms(), list)


if __name__ == '__main__':
    unittest.main()

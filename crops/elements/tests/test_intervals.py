"""Testing facilities for crops.elements.intervals"""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.elements import intervals as cei

import unittest

_INTEGER_1 = 5
_INTEGER_2 = 1
_INTEGER_3 = 20
_INTEGER_4 = 17
_INTEGER_LIST_1 = [5, 12]
_INTEGER_LIST_2 = [12, 5]
_INTEGER_LIST_3 = [1, 6]
_INTEGER_LIST_4 = [20, 22]
_INTEGER_LIST_5 = [11, 25]
_STRING_1 = 'test string'
_STRING_2 = 'test tag'

class TestCropsIntervals(unittest.TestCase):
    def test_intervalise_1(self):
        expected_subint = [[5, 5]]

        calc_subint = cei._intervalise(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intervalise_2(self):
        expected_subint = [[5, 5]]

        calc_interval = cei._intervalise(_INTEGER_1)
        calc_subint = cei._intervalise(calc_interval).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intervalise_3(self):
        expected_subint = [[5, 12]]

        calc_subint = cei._intervalise(_INTEGER_LIST_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intervalise_4(self):
        expected_subint = [[5, 12]]

        calc_subint = cei._intervalise(_INTEGER_LIST_2).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_union_1(self):
        expected_subint = [[5, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        calc_subint = interval.union(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_union_2(self):
        expected_subint = [[1, 1], [5, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        calc_subint = interval.union(_INTEGER_2).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_union_3(self):
        expected_subint = [[5, 12], [20, 2222220]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        calc_subint = interval.union(_INTEGER_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_union_4(self):
        expected_subint = [[1, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        calc_subint = interval.union(_INTEGER_LIST_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_union_5(self):
        expected_subint = [[5, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        calc_subint = interval.union(_INTEGER_LIST_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_1(self):
        expected_subint = [[5, 5]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_2(self):
        expected_subint = []

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_2).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_3(self):
        expected_subint = [[20, 20]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_4(self):
        expected_subint = []

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_5(self):
        expected_subint = [[5, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_LIST_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_6(self):
        expected_subint = [[5, 6]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_LIST_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_7(self):
        expected_subint = [[20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_LIST_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_intersection_8(self):
        expected_subint = [[11, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.intersection(_INTEGER_LIST_5).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_1(self):
        expected_subint = [[6, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_2(self):
        expected_subint = [[5, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_2).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_3(self):
        expected_subint = [[5, 12], [21, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_4(self):
        expected_subint = [[5, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_5(self):
        expected_subint = [[20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_LIST_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_6(self):
        expected_subint = [[7, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_LIST_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_7(self):
        expected_subint = [[5, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_LIST_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_subtract_8(self):
        expected_subint =[[5, 10]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.subtract(_INTEGER_LIST_5).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_1(self):
        expected_subint = [[6, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_2(self):
        expected_subint = [[1, 1], [5, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_2).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_3(self):
        expected_subint = [[5, 12], [21, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_4(self):
        expected_subint = [[5, 12], [17, 17], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_5(self):
        expected_subint = [[20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_LIST_1).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_6(self):
        expected_subint = [[1, 4], [7, 12], [20, 22]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_LIST_3).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_7(self):
        expected_subint = [[5, 12]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_LIST_4).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_symdiff_8(self):
        expected_subint = [[5, 10], [13, 19], [23, 25]]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.symdiff(_INTEGER_LIST_5).subint

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_terminals_1(self):
        expected_subint = [5, 12]

        interval = cei._intervalise(_INTEGER_LIST_1)

        calc_subint = interval.terminals()

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_terminals_2(self):
        expected_subint = [5, 22]

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_subint = interval.terminals()

        self.assertListEqual(calc_subint, expected_subint)

    def test_intinterval_n_elements_1(self):
        expected_number = 8

        interval = cei._intervalise(_INTEGER_LIST_1)

        calc_number = interval.n_elements()

        self.assertEqual(calc_number, expected_number)

    def test_intinterval_n_elements_2(self):
        expected_number = 11

        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        calc_number = interval.n_elements()

        self.assertEqual(calc_number, expected_number)

    def test_intinterval_contains_1(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        self.assertTrue(interval.contains(_INTEGER_1))

    def test_intinterval_contains_2(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        self.assertFalse(interval.contains(_INTEGER_2))

    def test_intinterval_contains_3(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        self.assertTrue(interval.contains(_INTEGER_LIST_1))

    def test_intinterval_contains_4(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        self.assertFalse(interval.contains(_INTEGER_LIST_3))

    def test_intinterval_contains_5(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        interval2 = cei._intervalise(_INTEGER_LIST_4)

        self.assertTrue(interval.contains(interval2))

    def test_intinterval_contains_6(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        interval2 = cei._intervalise(_INTEGER_LIST_5)

        self.assertFalse(interval.contains(interval2))

    def test_intinterval_description_1(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        self.assertTrue('description' in interval.tags)
        self.assertEqual(interval.tags['description'], 'intinterval')

    def test_intinterval_description_2(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        interval.description(_STRING_1)

        self.assertTrue('description' in interval.tags)
        self.assertEqual(interval.tags['description'], _STRING_1)

    def test_intinterval_addtag_1(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        interval.addtag(_STRING_2, _STRING_1)

        self.assertTrue(_STRING_2 in interval.tags)
        self.assertEqual(interval.tags[_STRING_2], _STRING_1)

    def test_intinterval_deltag_1(self):
        interval = cei._intervalise(_INTEGER_LIST_1)
        interval = interval.union(_INTEGER_LIST_4)

        interval.addtag(_STRING_2, _STRING_1)

        interval.deltag(_STRING_2)

        self.assertFalse(_STRING_2 in interval.tags)

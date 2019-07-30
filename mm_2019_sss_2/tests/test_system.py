"""
Unit and regression test for the mm_2019_sss_2 package.
"""

# Import package, test suite, and other packages as needed
import mm_2019_sss_2
import pytest
import sys

import numpy as np

def test_mm_2019_sss_2_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mm_2019_sss_2" in sys.modules

## Tests for system class
def test_initialize_random_simulation_np(box_length, num_particles):
        try:
            num_particles += 1
        except TypeError:
            print('Number of particles is not a numeric type (float or int).')

def test_initialize_random_simulation_bl(box_length, num_particles):
        try:
            box_length += 1
        except TypeError:
            print('Box length is not a numeric type (float or int).')

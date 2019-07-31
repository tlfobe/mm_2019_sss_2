"""
Unit tests for the energy.py methods
"""

import mm_2019_sss_2
import numpy as np
import pytest
import sys

def test_energy_factory_lj():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('LJ', epsilon = 0.5, sigma = 1)
    assert energy_function.calc_energy(2) == -0.03076171875

def test_energy_factory_buckingham():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('Buckingham', rho = 1, a = 1, c = 1)
    assert energy_function.calc_energy(1) == -0.6321205588285577

def test_energy_factory_unitlesslj():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('UnitlessLJ')
    assert energy_function.calc_energy(1) == 0.0

def test_energy_calculator():
    energy_obj = mm_2019_sss_2.energy.Energy()
    energy_obj.simulation_cutoff = 3
    assert(energy_obj.calculate_tail_correction(2, 10) == -0.0012405555234009783)
    test_coords = np.array([[0,0,0], [0,0,1.5]])
    assert(energy_obj.calculate_initial_energy(test_coords, 10) == -0.03059177374878155)
    energy_obj.simulation_cutoff = 1
    assert(energy_obj.calculate_initial_energy(test_coords, 10) == 0) 

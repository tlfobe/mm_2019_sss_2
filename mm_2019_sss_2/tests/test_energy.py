"""
Unit tests for the energy.py methods
"""

import mm_2019_sss_2
import pytest
import sys

def test_energy_factory_lj():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('LJ', epsilon = 0.5, sigma = 1)
    assert energy_function(2) == -0.03076171875

def test_energy_factory_buckingham():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('Buckingham', rho = 1, A = 1, C = 1)
    assert energy_function(1) == -0.6321205588285577

def test_energy_factory_unitlesslj():
    energy_factory = mm_2019_sss_2.energy.potentialEnergyFactory()
    energy_function = energy_factory.build_energy_method('UnitlessLJ')
    assert energy_function(1) == 0.0


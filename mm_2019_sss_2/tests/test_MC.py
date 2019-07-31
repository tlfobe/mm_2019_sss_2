"""
This is a unit test for monte_carlo.py.
"""

# Impoart package, test suit, and other packages as needed
import numpy as np
import mm_2019_sss_2 as MM2
import pytest
from argparse import Namespace

default_system = MM2.system.SystemSetup(num_particles=20, reduced_density=0.9)
MC_energy = MM2.energy.Energy()
default_namespace = Namespace(build_method='random', energy_function='UnitlessLJ', \
filename=None, freq=1000, max_displacement=0.1, n_steps=1000000, num_particles=20, \
reduced_density=0.9, reduced_temperature=0.9, simulation_cutoff=3.0, tune_displacement=True)

default_MC = MM2.monte_carlo.MonteCarlo(new_system=default_system, energy=MC_energy, arguments=default_namespace)
# Here we use default_system, MC_energy and default_namespace to generate an instance of the class MonteCarlo

@pytest.mark.parametrize("delta_e, beta, expected", [
    (-5.0, 1.0, True), (0.1, 1.0, True)
])

def test_accept_or_reject(delta_e, beta, expected):
    """
    For delta_e < 0, the expected result should always be true (proposed movement accepted).
    For delta_e > 0, the random seed 2019 makes the first random number always be 0.9034822144192743.
    In the test case, since we have delta_e = 0.1 and beta = 1.0, p_acc = 0.9048374180359595, which 
    is larger than the random number, so the expected result should be "True".
    """

    np.random.seed(2019)
    calculated = default_MC.accept_or_reject(delta_e, beta)
    try:
        assert expected == calculated
    finally:
        np.random.seed()

@pytest.mark.parametrize("n_trials, n_accept, max_displacement, expected", [
    (100, 36, 10, (8, 0, 0)), (100, 40, 10, (10, 0, 0)), (100, 44, 10, (12, 0, 0))
])

def test_adjust_displacement(n_trials, n_accept, max_displacement, expected):
    """
    First test case (100, 36, 10, (8, 0, 0)): acc_rate = 0.36 < 0.38, so max_displacement = 8
    Second test case (100, 40, 10, (10, 0, 0)): 0.38 < acc_rate = 0.40 < 0.42, so max_displacement = 10
    Third test case (100, 44, 10, (12, 0, 0)): acc_rate = 0.44 > 0.42, so max_displacement = 12
    """
    calculated = default_MC.adjust_displacement(max_displacement, n_accept, n_trials)
    assert expected == calculated

def test_parse_arguments_defaults():
    """
    This testing function tests if the defaults are overwritten correctly if the arugments are specified 
    differently from the defaults. Codewisely, note that default_MC.argument is a namespace, which can be 
    converted to a dictionary by using vars(default_MC.arguments). Then, by using list(vars(default_MC.arguments).values()), 
    we can get the list of defaults, which is to be compared with the list of specified values. 
    """

    # create an instance of MonteCarlo for the test case (all arugments are specified to be different from the defaults 
    # except for the build method)
    test_system = MM2.system.SystemSetup(num_particles=50, reduced_density=0.2)
    test_namespace = Namespace(build_method='file', energy_function='LJ', \
    filename=None, freq=500, max_displacement=0.5, n_steps=500000, num_particles=50, \
    reduced_density=0.2, reduced_temperature=0.5, simulation_cutoff=2.0, tune_displacement=False)
    test_energy = MM2.energy.Energy()
    test_case = MM2.monte_carlo.MonteCarlo(new_system=test_system, energy=test_energy, arguments=test_namespace)
    
    # things to be compared
    test_case_values = list(vars(test_case.arguments).values())
    default_values = list(vars(default_MC.arguments).values())
    
    # start to compare 
    count = 0
    for i in range(len(test_case_values)):
        print("test_case:", test_case_values[i])
        print("default_values:", default_values[i])
        if test_case_values[i] != default_values[i]:
            count += 1

    expected_count = len(test_case_values) - 1
    assert  expected_count == count







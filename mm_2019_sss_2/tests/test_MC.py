"""
This is a unit test for monte_carlo.py.
"""

# Impoart package, test suit, and other packages as needed
import numpy as np
import mm_2019_sss_2 as MM2

import pytest

MC_system = MM2.system.SystemSetup
MC_energy = MM2.energy.Energy
test_MC = MM2.monte_carlo.MonteCarlo(MC_system, MC_energy)

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
    calculated = test_MC.accept_or_reject(delta_e, beta)
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
    Second test case (100, 40, 10, (8, 0, 0)): 0.38 < acc_rate = 0.40 < 0.42, so max_displacement = 10
    Third test case (100, 44, 10, (8, 0, 0)): acc_rate = 0.44 > 0.42, so max_displacement = 12
    """
    calculated = my_MC.adjust_displacement(n_trials, n_accept, max_displacement)
    assert expected == calculated

def test_run_simulation




import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_initial_state(method='random', file_name=None, num_particles=None, \
 box_length=None):
    """Function generates initial coordinates for a LJ fluid simulation

    This function can read coordinates either from a file (NIST LJ Fluid \
     Configurations) or generates
    a random configuration.

    Parameters
    ----------
    method : str
        String the method to use to build the initial configuration for the LJ \
         fluid simulation. Possible values are 'random'  or 'file' (Default \
          value is 'random')
    file_name : str
        String of the the filename containing the initial starting coordinates. \
         Only required when using the 'fille' method (Default value = None)
    num_particles : int
        Number of particules to use when populating the simualtion box with the \
         'random' method (Default value = None)
    box_length : float
        Size of one vertices of the simulation box. (Default value = None)

    Returns
    -------
    coordinates : array_like
        A (num_particles x 3) numpy array containing the coordinates of each \
         LJ particle.

    Examples
    --------
    >>> generate_initial_state('random', num_particles = 1000, box_length = 20)
    array([[ 1.10202674,  4.24975121, -5.03322129],
        [ 9.13676284,  4.78807621, -8.26008762],
        [ 6.24720765, -7.17769567,  9.61620896],
        ...,
        [-3.47864571,  2.32867699, -1.31176807],
        [ 1.3302019 , -3.4160087 , -1.34698966],
        [ 0.56410479, -1.2309513 ,  4.71009776]])
    """

    if method is 'random':

        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length

    elif method is 'file':

        coordinates = np.loadtxt(file_name, skiprows=2, usecols=(1,2,3))

    return coordinates


def lennard_jones_potential(rij2):
    """Function evaluates the unitless LJ potential given a squared distance

    Parameters
    ----------
    rij2 : float
        Distance squared between two particles

    Returns
    -------
    value : float
        Unitless LJ potential energy
    """

    sig_by_r6 = np.power(1 / rij2, 3)
    sig_by_r12 = np.power(sig_by_r6, 2)

    return 4.0 * (sig_by_r12 - sig_by_r6)

def calculate_tail_correction(box_length, cutoff, number_particles):
    """This function computes the standard tail energy correction for the LJ \
     potential

    Parameters
    ----------
    box_length : float/int
        length of simulation box
    cutoff: float/int
        the cutoff for the tail energy truncation
    num_particles: int
        number of particles

    Returns
    ------_
    e_correction: float
        tail correction of energy
    """

    volume = np.power(box_length, 3)
    sig_by_cutoff3 = np.power(1.0 / cutoff, 3)
    sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
    e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3

    e_correction *= 8.0 / 9.0 * np.pi * number_particles / volume * \
     number_particles

    return e_correction

def minimum_image_distance(r_i, r_j, box_length):
    # This function computes the minimum image distance between two particles

    rij = r_i - r_j
    rij = rij - box_length * np.round(rij / box_length)
    rij2 = np.dot(rij, rij)
    return rij2

def get_particle_energy(coordinates, box_length, i_particle, cutoff2):

    """This function computes the minimum image distance between two particles

    Parameters
    ----------
    r_i: list/array
        the potitional vection of the particle i
    r_j: list/array
        the potitional vection of the particle j
    box_length : float/int
        length of simulation box

    Returns
    ------_
    rij2: float
        the square of the shortest distance between the two particles and their \
         images
    """

    e_total = 0.0

    i_position = coordinates[i_particle]

    particle_count = len(coordinates)

    for j_particle in range(particle_count):

        if i_particle != j_particle:

            j_position = coordinates[j_particle]

            rij2 = minimum_image_distance(i_position, j_position, box_length)

            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total


def calculate_total_pair_energy(coordinates, box_length, cutoff2):
    """This function computes the sum of all pairwise VDW energy between each \
     pair of particles in the system.

    Parameters
    ----------
    coordinates : np.array
        An array of atomic coordinates. Size should be (n, 3) where n is the \
         number of particles.
    box_length : float
        A float indicating the size of the simulation box. Can be either \
         hard-coded or calculated using num_particles and reduced_density.
    cutoff2: float
        The square of the simulation_cutoff, which is the cutoff distance \
         between two interacting particles.

    Returns
    -------
    e_total : float
        The sum of all pairwise VDW energy between each pair of particles in \
         the system.
    """

    e_total = 0.0
    particle_count = len(coordinates)

    for i_particle in range(particle_count):
        for j_particle in range(i_particle):

            r_i = coordinates[i_particle]
            r_j = coordinates[j_particle]
            rij2 = minimum_image_distance(r_i, r_j, box_length)
            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total

def accept_or_reject(delta_e, beta):
    """Accept or reject a move based on the energy difference and system \
     temperature.

    This function uses a random numbers to adjust the acceptance criteria.

     Parameters
     ----------
     delta_e : float
        The difference between the proposed and current energies.
     beta : float
        The inverse value of the reduced temperature.

     Returns
     -------
     accept : booleen
        Either a "True" or "False" to determine whether to reject the trial.
    """

    if delta_e < 0.0:
        accept = True

    else:
        random_number = np.random.rand(1)
        p_acc = np.exp(-beta * delta_e)

        if random_number < p_acc:
            accept = True
        else:
            accept = False

    return accept


def adjust_displacement(n_trials, n_accept, max_displacement):
    """Change the acceptance criteria to get the desired rate.

    When the acceptance rate is too high, the maximum displacement is adjusted \
     to be higher.
    When the acceptance rate is too low, the maximum displacement is \
     adjusted lower.

    Parameters
    ----------
    n_trials : integer
        The number of trials that have been performed when the function is \
         initiated.
    n_accept : integer
        The current number of accepted trials when the function is initiated.
    max_displacement : float
        The specified maximum value for the displacement of the trial.

    Returns
    -------
    max_displacement : float
        The adjusted displacement based on the acceptance rate.
    n_trials : integer, 0
        The new number of trials.
    n_accept : integer, 0
        The new number of trials.
    """
    acc_rate = float(n_accept) / float(n_trials)
    if (acc_rate < 0.38):
        max_displacement *= 0.8

    elif (acc_rate > 0.42):
        max_displacement *= 1.2

    n_trials = 0
    n_accept = 0

    return max_displacement, n_trials, n_accept

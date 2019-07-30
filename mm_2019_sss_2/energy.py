import numpy as np
from abc import ABC, abstractmethod


class EnergyFunction(ABC):
    @abstractmethod
    def calc_energy(self):
        pass

    def cutoff_correction(self):
        pass


class LJ(EnergyFunction):
    """
    Set-up for the Lennard--Jones potential.

    Parameters
    ----------
    epsilon: float, int

    sigma: float, int


    Returns
    --------
    A class?
    """
    def __init__(self, epsilon, sigma):
        self.sigma = sigma
        self.epsilon = epsilon

    def calc_energy(self, r):
        return (4 * self.epsilon * ((self.sigma / r) ** 12
                                    - (self.sigma / r) ** 6))

    def cutoff_correction(self, cutoff, number_particles, box_length):
        pass


class Buckingham(EnergyFunction):
    """
    Set-up for the Buckingham potential.

    Parameters
    ----------
    epsilon: float, int

    sigma: float, int


    Returns
    --------
    A class?
    """
    def __init__(self, rho, A, C):
        self.rho = rho
        self.A = A
        self.C = C

    def calc_energy(self, r):
        return (self.A * np.exp(-r / self.rho) - self.C / r ** 6)

    def cutoff_correction(self, cutoff, number_particles, box_length):
        pass


class UnitlessLJ(EnergyFunction):
    def __init__(self):
        pass

    def calc_energy(self, r):
        return (4.0 * (np.power(1 / r, 12) - np.power(1 / r, 6)))

    def cutoff_correction(self, cutoff, number_particles, box_length):
        # This function computes the standard tail energy correction for the LJ potential

        volume = np.power(box_length, 3)
        sig_by_cutoff3 = np.power(1.0 / cutoff, 3)
        sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
        e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3

        e_correction *= 8.0 / 9.0 * np.pi * number_particles / volume * number_particles

        return e_correction


class potentialEnergyFactory:
    def __init__(self):
        self.methods = {'LJ': LJ,
                        'Buckingham': Buckingham,
                        'UnitlessLJ': UnitlessLJ,
                        }

    def build_energy_method(self, potential_type, **kwargs):
        energy_class = self.methods[potential_type](**kwargs)
        return (energy_class)


class Energy:
    def __init__(self, potential_type='UnitlessLJ', simulation_cutoff=3.0,
                 **kwargs):
        self.energy_obj = potentialEnergyFactory().build_energy_method(
            potential_type, **kwargs)
        self.simulation_cutoff = simulation_cutoff

    def calculate_tail_correction(self, number_particles, box_length):
        """
        This function computes the standard tail energy correction for the LJ potential

        Parameters
        ----------
        box_length : float, int
            length of of side of the simulation box (cube)
        cutoff: float, int
            the cutoff for the tail energy truncation
        num_particles: int
            number of particles

        Return
        ------
        e_correction: float
            tail correction of energy
        """
        e_correction = self.energy_obj.cutoff_correction(
            self.simulation_cutoff, number_particles, box_length)
        return e_correction

    def _minimum_image_distance(self, r_i, r_j, box_length):
        """
        Calculates the shortest distance between a particle and another
        instance in a periodic boundary condition image

        Parameters
        ----------
        r_i: np.array([n,3])
            The x, y, z coordinates for a particle, i.
        r_j: np.array([n,3])
            The x, y, z coordinates for a particle, j.
        box_length: float, int
            The length of a side of the side box for the periodic boundary.

        Returns
        -------
        rij2: np.array([n,3])
            The minimum image distance between the two particles, r_i and r_j.
        """
        # This function computes the minimum image distance between two particles
        rij = r_i - r_j
        rij = rij - box_length * np.round(rij / box_length)
        rij2 = np.dot(rij, rij)
        return rij2

    def calculate_initial_energy(self, coordinates, box_length):
        """Iterates over a set of coordinates to calculate total system energy

        This function computes the sum of all pairwise VDW energy between each
        pair of particles in the system. This is the first instance of the
        energy calculation. Subsequent uses call calculate_pair_energy.

        Parameters
        ----------
        coordinates : np.array([n,3])
            An array of atomic coordinates. Size should be [n, 3] where n is the
            number of particles.
        box_length : float
            A float indicating the size of the simulation box. Can be either
            hard-coded or calculated using num_particles and reduced_density.
        cutoff: float
            The square of the simulation_cutoff, which is the cutoff distance
            between two interacting particles.

        i_particle: int
            Intitial particle for pairwise count

        Returns
        -------
        e_total : float
            The sum of all pairwise VDW energy between each pair of particles in
            the system.

        e_total = 0.0
        particle_count = len(coordinates)
        for i_particle in range(particle_count):
            for j_particle in range(i_particle):
                r_i = coordinates[i_particle]
                r_j = coordinates[j_particle]
                rij2 = self._minimum_image_distance(r_i, r_j, box_length)
                print(f'rij2: {rij2}')
                if rij2 < self.simulation_cutoff ** 2:
                    e_pair = self.energy_obj.calc_energy(np.sqrt(rij2))
                    print(f'e pair: {e_pair}')
                    e_total += e_pair
        return e_total


    def calculate_pair_energy(self, coordinates, box_length, i_particle):
    """This function computes the sum of all pairwise VDW energy between each
        pair of particles in the system.

    Parameters
    ----------
    coordinates : np.array
        An array of atomic coordinates. Size should be (n, 3) where n is the
        number of particles.
    box_length : float
        A float indicating the size of the simulation box. Can be either
        hard-coded or calculated using num_particles and reduced_density.
    cutoff: float
        The square of the simulation_cutoff, which is the cutoff distance between
        two interacting particles.

    i_particle: integer
        Intitial particle for pairwise count

    Returns
    -------
    e_total : float
        The sum of all pairwise VDW energy between each pair of particles in
        the system.
    """

        # This function computes the energy of a particle with
        # the rest of the system

        e_total = 0.0

        i_position = coordinates[i_particle]

        particle_count = len(coordinates)

        for j_particle in range(particle_count):

            if i_particle != j_particle:

                j_position = coordinates[j_particle]

                rij2 = self._minimum_image_distance(i_position, j_position,
                                                    box_length)

                if rij2 < self.simulation_cutoff ** 2:
                    e_pair = self.energy_obj.calc_energy(np.sqrt(rij2))
                    e_total += e_pair
        return e_total


def main():
    energy_factory = potentialEnergyFactory()
    lj_energy = energy_factory.build_energy_method('UnitlessLJ')
    print(lj_energy.calc_energy(2.0))


if __name__ == "__main__":
    main()

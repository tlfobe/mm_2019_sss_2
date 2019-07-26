import numpy as np
from abc import ABC, abstractmethod

class EnergyFunction(ABC):
    @abstractmethod
    def calc_energy(self):
        pass
    def cutoff_correction(self):
        pass


class LJ:
    def __init__(self, epsilon, sigma):
        self.sigma = sigma
        self.epsilon = epsilon
    
    def calc_energy(self, r):
        return(4 * self.epsilon * ( (self.sigma / r)**12 - (self.sigma /     r)**6 ))

    def cutoff_correction(self, cutoff, number_particles, box_length):
        pass

class Buckingham:
    def __init__(self, rho, A, C):
        self.rho = rho
        self.A = A
        self.C = C

    def calc_energy(self, r):
        return( self.A * np.exp( -r / self.rho) - self.C / r**6 )

    def cutoff_correction(self, cutoff, number_particles, box_length):
        pass

class UnitlessLJ:
    def __init__(self):
        pass
    
    def calc_energy(self, r):
        return(4.0 * (np.power(1/r , 6) - np.power(1/r, 12)))

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
        self.methods = {'LJ'        : LJ,
                        'Buckingham': Buckingham,
                        'UnitlessLJ': UnitlessLJ,
                        }

    def build_energy_method(self, potential_type, **kwargs):
        energy_class = self.methods[potential_type](**kwargs)
        return (energy_class.calc_energy)


class EnergyCalculator:
    def __init__(self, potential_type = 'UnitlessLJ', **kwargs):
        self.energy_function = potentialEnergyFactory().build_energy_method(potential_type, **kwargs)


    def _minimum_image_distance(r_i, r_j, box_length):
        # This function computes the minimum image distance between two particles
    ​
        rij = r_i - r_j
        rij = rij - box_length * np.round(rij / box_length)
        rij2 = np.dot(rij, rij)
        return rij2



    def calculate_initial_energy(self, coordinates, box_length, cutoff):
        e_total = 0.0
        e_total += self.energy_function.cutoff_correction(cutoff, len(coordinates), box_length)
        particle_count = len(coordinates)
        for i_particle in range(particle_count):
            for j_particle in range(i_particle):
                r_i = coordinates[i_particle]
                r_j = coordinates[j_particle]
                rij2 = _minimum_image_distance(r_i, r_j, box_length)
                if rij2 < cutoff**2:
                    e_pair = self.energy_function(rij2)
                    e_total += e_pair
        return e_total


    def calculate_pair_energy(self, coordinates, box_length, cutoff, i_particle):

        #This function computes the energy of a particle with
        #the rest of the system

        e_total = 0.0

        i_position = coordinates[i_particle]

        particle_count = len(coordinates)

        for j_particle in range(particle_count):

            if i_particle != j_particle:
                
                j_position = coordinates[j_particle]
    
                rij2 = minimum_image_distance(i_position, j_position, box_length)
    
                if rij2 < cutoff**2:
                    e_pair = self.energy_function(rij2)
                    e_total += e_pair
        return e_total

    



def main():
    energy_factory = potentialEnergyFactory()
    lj_energy = energy_factory.build_energy_method('UnitlessLJ')
    print(lj_energy(2.0))

if __name__ == "__main__":
    main()
    




import numpy as np
from abc import ABC, abstractmethod

class EnergyFunction(ABC):
    @abstractmethod
    def calc_energy(self):
        pass


class LJ:
    def __init__(self, epsilon, sigma):
        self.sigma = sigma
        self.epsilon = epsilon
    
    def calc_energy(self, rij):
        return(4 * self.epsilon * ( (self.sigma / r)**12 - (self.sigma /     r)**6 ))

class Buckingham:
    def __init__(self, rho, A, C):
        self.rho = rho
        self.A = A
        self.C = C

    def calc_energy(self, r):
        return( self.A * np.exp( -r / self.rho) - self.C / r**6 )

class UnitlessLJ:
    def __init__(self):
        pass
    
    def calc_energy(self, r):
        return(4.0 * (np.power(1/r , 6) - np.power(1/r, 12)))

class potentialEnergyFactory:
    def __init__(self):
        self.methods = {'LJ'        : LJ,
                        'Buckingham': Buckingham,
                        'UnitlessLJ': UnitlessLJ,
                        }

    def build_energy_method(self, potential_type, **kwargs):
        energy_class = self.methods[potential_type](**kwargs)
        return (energy_class.calc_energy)

def main():
    energy_factory = potentialEnergyFactory()
    lj_energy = energy_factory.build_energy_method('UnitlessLJ')
    print(lj_energy(2.0))

if __name__ == "__main__":
    main()
    




import numpy as np

class SystemSetup:
    """A class object that initializes the system for your Monte Carlo
    calculation.

    By defining the method as either 'random' or 'file' you an either generate
    a random box of Ar atoms, or you can read in the coordinates the specifying
    the file using the 'filename' keyword.

    Parameters
    ----------
    method : str
        The method is either 'random' or 'file'. By default the method is set
        to 'random'.
    num_particles : int
        The number of particles should be defined as an integer value. This
        keyword is only really needs to be defined if you specify the method
        to be 'random'. By default the value is set to 20.
    box_length : float
        The box is assumed to be cubic. This keyword is also only needed when
        you specify the method to be 'random'. By default the value is set
        to 3.0 Angstroms.
    filename : str
        The file name if the method is set to file. By default the system is
        initialized using random coordinates. If you want to read in a file,
        specify the method to be 'file' and provide a filename as a string
        using this keyword.

    Scenario 1: Generating a random box with coordinates.
    >>> test_system_1 = SystemSetup(method="random")
    >>> print(test_system_1.n_particles)
    >>> print(test_system_1.coordinates)

    Scenario 2: Read a file that initializes your system.
    >>> test_system_2 = SystemSetup(method="file", filename="input.dat")
    >>> print(test_system_2.n_particles)
    >>> print(test_system_2.coordinates)
    """

    def __init__(self, method: str = 'random', num_particles: int = 20,
                 box_length: (int, float) = 3.0, filename: str = None):

        self.method = method
        self.box_length = float(box_length)

        if self.method == 'random':
            self._initialize_random_simulation_(box_length, num_particles)
        elif self.method == 'file':
            try:
                self._read_info_from_file_(filename)
            except FileNotFoundError:
                filename = input("Couldn't find your filename, try again: ")
                self._read_info_from_file_(filename)
                try:
                    self._read_info_from_file_(filename)
                except FileNotFoundError:
                    self._read_in_error_(num_particles, box_length)
        else:
            raise TypeError('You are using a method that is not supported '
                            'at this moment.')

    def _read_in_error_(self, num_particles, box_length):
        print('Either you entered the incorrect file or the file was not '
              'found.')
        print('Initializing a Monte Carlo simulation with the following '
              'parameters...')
        print(f'Number of particles: {num_particles}')
        print(f'Box length: {box_length} Angstroms')
        self._initialize_random_simulation_(box_length, num_particles)

    def _read_info_from_file_(self, filename):
        self.coordinates = np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3))
        self.n_particles = len(self.coordinates)

    def _initialize_random_simulation_(self, box_length, num_particles):
        self.n_particles = num_particles
        self.coordinates = (0.5 - np.random.rand(num_particles, 3)) * \
                           box_length

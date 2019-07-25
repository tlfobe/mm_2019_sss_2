import numpy as np
import os

class System:
	"""This is a class object that initializes a system for your Monte Carlo calculation.

	Parameters
	----------
	method : str
		The method is either 'random' or 'file'. By default the method is set to 'random'.
	num_particles : int
		The number of particles should be defined as an integer value. By default the value is 20.
	box_length : float
		The box is assumed to be cubic. By default the value is set to 3.0 Angstroms.
	filename : str
		The file name if the method is set to file.
	"""

	def __init__(self, method='random', num_particles=20, box_length=3.0, filename=None):
		self.box_length = box_length
		self.method = method

		if self.method == 'random':
			self._initialize_random_simulation_(box_length, num_particles)

		elif self.method == 'file':
			try:
				self._read_info_from_file_(filename)
			except FileNotFoundError:
				self._read_in_error(num_particles, box_length)
		else:
			raise TypeError('You are using a method that is not supported at  this moment.')
	
	def _read_in_error(self, num_particles, box_length):
		print('Initializing a Monte Carlo simulation with: ')
		print('Either you entered the incorrect file or the file was not found.')
		print(f'Number of particles: {num_particles}')
		print(f'Box length: {box_length} Angstroms')
		self._initialize_random_simulation_(box_length, num_particles)

	def _read_info_from_file_(self, filename):
		self.coordinates = np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3))
		self.n_particles = len(self.coordinates)

	def _initialize_random_simulation_(self, box_length, num_particles):
		self.n_particles = num_particles
		self.coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length

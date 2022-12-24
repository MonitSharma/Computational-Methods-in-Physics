'''
Suppose we want to simulate the motion and dynamic interactions between n particles in a closed 2-dimensional container.
The task at hand is to numerically integrate Newton’s equations for a set of particles, with the interaction force between
any two particles being derived from the Lennard Jones potential. The net force on each particle is then to be used to evaluate
the system over time, using the Velocity Verlet algorithm for a timestep Delta t.

The program was essentially implemented using two functions, one to calculate the resultant force on each particle and the other to
carry out the Velocity Verlet algorithm. Heavy use of Numpy arrays and advanced indexing was made, effectively resulting in an implementation
without the use of any looping constructs (besides for plotting). In order to reduce the number of calculations needed to find all of
the forces, an upper triangular matrix of relative position vectors was constructed, and the force was then calculated on the position 
vectors within this upper triangular matrix.


Using Newton’s Third Law of Motion, the corresponding lower triangular matrix could be found (by using the results in the upper
triangular matrix rather than carrying out another set of calculations). Summing these two matrices resulted in a matrix with the
component forces for the i th particle in the i th row. Thus by summing over the 1st axis, the resultant force for the i th particle
was found, where i is in {1, ..., n}. To further reduce the number of calculations, a norm cutoff parameter was added.

A number of readily available Numpy functions were used for efficiency, including: np.transpose for transposing arrays, np.linalg.norm
for finding the position vectors, np.random.uniform for generating the particle positions and velocities, as well as the Numpy functions
for array addition and multiplication. Indexing arrays were also generated using np.greater and np.less.


The same holds for the Velocity Velocity algorithm - the required equations were essentially implemented as the addition of arrays,
eliminating the need for any for loops. Within the Velocity Verlet algorithm a set of boundary reflection conditions were implemented 
for particles crossing the container boundary - effectively the particles were reflected off of the boundary, with the assumption of 
conservation of energy.

Lastly, for the graphical aspect of the program the scatter functionality provided by matplotlib was used to plot the position vectors 
for each particle with each iteration of the while loop.
'''



import numpy as np
from matplotlib import pyplot as plt

# Constants
length = 50
v_max = 15
delta_t = 0.01
r_cut = 10
particle_count = 9

# Initialise simulation
positions = np.random.uniform(0, length, (particle_count, 2))
velocities = np.random.uniform(-v_max, v_max, (particle_count, 2))
forces = np.zeros((particle_count, 2), dtype=np.float64)

# Function to calculate the forces using the Lennard-Jones Potential
def forces_calc():

    # Create square matrices of x,y coords for particles for i > j
    pos_tiled_x = np.tile(positions[:, 0], (particle_count, 1))
    pos_tiled_y = np.tile(positions[:, 1], (particle_count, 1))

    r_ij = np.empty((particle_count, particle_count, 2), dtype=np.float64)
    # Finding the position vectors r_ij
    r_ij[:, :, 0] = np.subtract(np.triu(np.transpose(pos_tiled_x)), np.triu(pos_tiled_x))
    r_ij[:, :, 1] = np.subtract(np.triu(np.transpose(pos_tiled_y)), np.triu(pos_tiled_y))
    norm_r_ij = np.linalg.norm(r_ij, axis=2)  # Finding the norm of the position vectors

    index_array = np.less(norm_r_ij, r_cut)  # Indexing array for particles satisfying the threshold condition
    calc_forces = np.zeros((particle_count, particle_count), dtype=np.float64)
    calc_forces[index_array] = np.subtract(np.multiply(np.power(norm_r_ij[index_array], -14), 48),
                                           np.multiply(np.power(norm_r_ij[index_array], -8), 24))  # Potential for ij
    calc_forces[np.isnan(calc_forces)] = 0.0  # Handling NaNs in the resulting upper triangular matrix
    r_ij[:, :, 0] = np.multiply(r_ij[:, :, 0], calc_forces)  # Resulting force in the x component for ij
    r_ij[:, :, 1] = np.multiply(r_ij[:, :, 1], calc_forces)  # Resulting force in the y component for ij
    r_ij[:, :, 0] = np.add(r_ij[:, :, 0], -np.transpose(r_ij[:, :, 0]))  # Corresponding lower triangular entry for x
    r_ij[:, :, 1] = np.add(r_ij[:, :, 1], -np.transpose(r_ij[:, :, 1]))  # Corresponding lower triangular entry for y

    global forces
    forces = np.sum(r_ij, axis=1)  # Resultant net force on each particle i

# Function to evaluate the motion of the particles
def velocity_verlet():

    global positions, velocities
    positions = np.add(np.add(positions, np.multiply(delta_t, velocities)), np.multiply(0.5*(delta_t**2), forces))  # Update position of particles

    # Boundary reflection conditions
    index_array_L = np.greater(positions, length)
    positions[index_array_L] = np.add(-positions[:, 0:2][index_array_L], 2*length)  # For beyond length L
    index_array_0 = np.less(positions, 0)
    positions[index_array_0] = -positions[index_array_0]  # For less than 0

    forces_prev = np.copy(forces)  # Copy forces from previous iteration
    forces_calc()  # Calculate the new forces

    velocities = np.add(velocities, np.multiply(delta_t/2, np.add(forces_prev, forces)))  # Update velocities

    # Boundary velocity inversion conditions
    velocities[np.logical_or(index_array_L, index_array_0)] = -velocities[np.logical_or(index_array_L, index_array_0)]

# Main code
plt.ion()
plt.show()

forces_calc()  # Calculating the initial forces to start the simulation

while True:
    velocity_verlet()
    plt.xlim((0, length))
    plt.ylim((0, length))
    plt.autoscale(False)
    plt.scatter(positions[:, 0], positions[:, 1])  # Plotting
    plt.pause(0.001)
    plt.clf()
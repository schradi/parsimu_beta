/*
 * SimProcess.h
 *
 *  Created on: 21.04.2015
 *      Author: jonas
 */

#ifndef SIMPROCESS_H_
#define SIMPROCESS_H_

#include "cell.h"
#include "particle.h"
#include "defines.h"
#include <math.h>
#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>


/**
 * Class Description SimProcess
 * todo
 */
class SimProcess{
public:
	real cell_size[DIM];			/**< +Width of each cell*/
	int  global_nc[DIM];			/**< +Global number of cells*/
	real global_size[DIM];			/**< +Size of the simulated area*/
	int  ic_start[DIM];				/**< +lowest cell-index in this process*/
	int  ic_stop[DIM];				/**< +highest cell-index in this process*/
	int  ip[DIM];					/**< +Position of this Process in map of Processes*/
	int  local_nc[DIM];				/**< +Number of cells within this Process*/
	real local_size[DIM];			/**< +Size of the simulated area in this Process*/
	int  neigh_lower[DIM];			/**< +Rank of the lower neighboring processes in each direction*/
	int	 neigh_upper[DIM];			/**< +Rank of the lower neighboring processes in each direction*/
	int  global_np[DIM];			/**< +Number of running Processes*/
	int  num_part;					/**< +Number of particles within this process*/
	int	 rank;						/**< +Rank of this Process*/
	real r_cut;						/**< +cutting-radius for short-range forces*/
	real start[DIM];				/**< +Starting-Point of this Process (left down point)*/
	real max_V;						/**< maximum of speed the particles are allowed to achieve*/
	int output_resolution;
	real t;

	/**
	 * Communication between the processes
	 */
	void communicate(Cell* cells);

	/**
	 * Communication between the processes in a specific Dimension
	 * @param com_d specific dimension
	 */
	void communicate(int com_d, Cell* cells);

	/**
	 * Calculates the Force between all Particles within the Cells of this Process
	 */
	void compF(Cell* cells);

	/**
	 * Calculates the Force between two specific Particles
	 * The force will be added to the force saved in the particle
	 * @param p_i Particle no 1
	 * @param p_j Particle no 2
	 */
	void compF(Particle* p_i, Particle* p_j);

	/**
	 * Calculates the new Velocity of the particles caused by the new force
	 * @param delta_t is the lenght of the timestep
	 */
	void compV(real delta_t, Cell* cells);

	/**
	 * Calculates the new Position of the particle caused by the new velocity
	 * @param delta_t is the lenght of the timestep
	 */
	void compX(real delta_t, Cell* cells);

	/**
	 * Calculates the force caused by the Lennard-Jones-Potential in a distance of p_X
	 * @param p_X Position of which the responsible process in needed
	 * @return rank of the process responsible for the Position p_X
	 */
	void force(real* p_X);

	/**
	 * Returns the id of the cell responsible for the Position p_X
	 * @param p_X Position of which the responsible process in needed
	 * @param id where the new Position will be stored
	 */
	int getLocalCellId(real* p_X);

	/**
	 * Create a new Particle
	 * @param p place, where new particle will be stored
	 */


	/**
	 * Calculates the global index out of a global position
	 * @param p current position
	 * @param n	global size
	 */
	int  index(int* p, int* n);

	/**
	 * Include Data
	 */
	void initData (Cell* cells); // initialize Data Structure and include Data

	/**
	 * Calculates the local index out of a global position
	 * @param p current position
	 * @param n	global size
	 */
	int  local_index(int* p, int* n);

	void matchOutputs();

	/**
	 * Puts particles in right cells
	 */
	void moveParticles(Cell* cells);

	/**
	 * Collects all outputs of the single processes in one process
	 * Process 0 creates the final output
	 */
	void gather_outputs();

	/**
	 * Saves all needed particleinformation calculated within this process
	 */
	void output(real delta_t, Cell* cells);

	void create_cells(Cell* p_cells);

	/**
	 * Constructor
	 */
	SimProcess (real p_r_cut, real* p_global_size, int* p_global_np);

	/**
	 * Destructor
	 */
	~SimProcess();

	/**
	 * Prints a status during the calculation
	 * @param c output textline
	 */
	void status_output(std::string c);

	/**
	 * Calculates the next timestep
	*/
	void timeIntegration(real delta_t, real t_end, Cell* cells);

	/**
	 * All Particles stored in adding are put in pl, adding = NULL
	 */
	void adding_to_pl(Cell* cells);

	/**
	 * Collect all pl of cells in a range and convert them to an array of double
	 */
	void code_pl(real* send_pl, int* icr_start, int* icr_stop, int com_d, int oth_d, int size, Cell* cells);

	/**
	* Convert an array of double into Particles and sort them in cells in a specified range
	*/
	void uncodePl(real* recv_pl, int* icr_start, int*icr_stop, int com_d, int oth_d, int length_recv, int size, Cell* cells);

	/**
	 * Send and receive a converted particle Cell to all neighboring processes in a direction
	 * @param com_d specific direction
	 * @param send_pl converted List of particles to be send
	 * @param recv_pl empty converted List of particles where the received particles will be stored
	 */
	void pl_send_recv(real* send_pl, int send_pl_l, real* recv_pl, int recv_pl_l, int com_d);

	int get_num_p(int* icr_start, int* icr_stop, int com_d, int oth_d, Cell* cells);

	void pl_l_send_recv(int* send_pl_l, int* recv_pl_l, int com_d);
};

#endif /* SIMPROCESS_H_ */

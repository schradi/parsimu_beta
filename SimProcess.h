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
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "timers.h"


/**
 * Class Description SimProcess
 * todo
 */
class SimProcess{
public:
	real lj_force_r_cut;
	real cell_size[DIM];			/**< +Width of each cell*/
	real delta_t;
	real t_end;
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
	long int global_num_part;
	long int  num_part;	/**< +Number of particles within this process*/
	long int	 num_ghost_part;
	int	 rank;						/**< +Rank of this Process*/
	real r_cut;						/**< +cutting-radius for short-range forces*/
	real start[DIM];				/**< +Starting-Point of this Process (left down point)*/
	real max_V;				/**< +maximum of speed the particles are allowed to achieve*/
	real sigma;					/**< Parameter sigma for the Lennard-Jones-Potential*/
	real epsilon; 				/**< Parameter epsilon for the Lennard-Jones-Potential*/
	TimerList* timer_list;
	real Vvar;
	int output_resolution;
	real t;
	int np;
	char* output_folder;
	bool log_time;
	bool log_energy;
	bool log_positions;
	bool log_velocity;
	bool log_id;
	bool folder_output;
	bool FALSCH;
	real E_kin;
	real E_pot;

	void copy_border_cells(int com_d, Cell* cells);

	void create_particles(ParticleList* new_pl);
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
	void compA(Cell* cells);

	/**
	 * Calculates the Force between two specific Particles
	 * The force will be added to the force saved in the particle
	 * @param p_i Particle no 1
	 * @param p_j Particle no 2
	 */
	void compA(Particle* p_i, Particle* p_j);

	/**
	 * Calculates the Force between all Particles within the Cells of this Process
	 */
	void compE(Cell* cells);

	/**
	 * Calculates the new Velocity of the particles caused by the new force
	 * @param delta_t is the lenght of the timestep
	 */
	void compV(Cell* cells);

	/**
	 * Calculates the new Position of the particle caused by the new velocity
	 * @param delta_t is the lenght of the timestep
	 */
	void compX(Cell* cells);

	/**
	 * Calculates the force caused by the Lennard-Jones-Potential in a distance of p_X
	 * @param p_X Position of which the responsible process in needed
	 * @return rank of the process responsible for the Position p_X
	 */
	void force(real* p_X, real* F);

	real lj_force(real r);

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
	void initData (Cell* cells, real* p_map); // initialize Data Structure and include Data

	void devide_symetric(real* p_map);
	void spread_global_info();
	void spread_local_info(real* p_map);
	void calculate_local_constants();

	/**
	 * Calculates the local index out of a global position
	 * @param p current position
	 * @param n	global size
	 */
	int  local_index(int* p);

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
	void output(Cell* cells, int outp_nr);

	void create_cells(Cell* p_cells);

	/**
	 * Constructor
	 */
	SimProcess(char* p_output_folder);

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
	void timeIntegration(Cell* cells, real* p_map);

	/**
	 * All Particles stored in adding are put in pl, adding = NULL
	 */
	void adding_to_pl(Cell* cells);

	/**
	 * Collect all pl of cells in a range and convert them to an array of double
	 */
	void code_range(real* send_pl, int* icr_start, int* icr_stop, Cell* cells);

	void code_p(Particle* p, real* code);

	void sort_particles(ParticleList* pl, real* p_map, ParticleList** sorted_particles, long int** pl_l);

	/**
	* Convert an array of double into Particles and sort them in cells in a specified range
	*/
	void uncode_in_range(real* recv_pl, int* icr_start, int*icr_stop, long int length_recv, Cell* cells, int pb_corr, int com_d);

	void uncode_in_range(real* recv_pl, int* icr_start, int*icr_stop, long int length_recv, Cell* cells){
		uncode_in_range(recv_pl, icr_start, icr_stop, length_recv, cells, 0, -1);
	};

	void uncode_p(real* code, Particle* p);

	void delete_pl(int* icr_start, int* icr_stop, Cell* cells);

	int get_num_p(int* icr_start, int* icr_stop, Cell* cells);

	void spread_particles(ParticleList* pl, Cell* cells, real* p_map);

	void create_particle_cloud(ParticleList* new_pl, real* r_start, real* r_stop, real resolution, real* p_V);
};

#endif /* SIMPROCESS_H_ */

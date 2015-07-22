/*
 * cell.h
 *
 *  Created on: 24.03.2015
 *      Author: jonas
 */

#ifndef CELL_H_
#define CELL_H_

#include"defines.h"
#include"particle.h"

class Cell{
public:
	ParticleList* adding; 	/**< If Cell is part of the outer border this contains all Particles that entered the cell */
	real cell_size[DIM];	/**< Cell-Width */
	int id;					/**< cell-id*/
	int num_part;			/**< */
	ParticleList* pl;		/**< List of Particles in the Cell */
	real start[DIM];		/**< Starting-Point of the cell */

	Cell(real* p_cell_size, int p_id, real* p_start);

	Cell(){num_part=0; adding=NULL; pl=NULL;};

	/*
	* Move all particles stored in adding to pl
	 * Adding is empty afterwards
	 */
	void adding_to_pl();

	/**
	 * Copy a existing Particle in the cell
	 * @param p_p pointer to existin particle
	 */
	void insertParticle(Particle* p_p);

	void insertParticle(ParticleList* p_pl);

	/**
	 * deletes the whole ParicleList pl
	 */
	void deletePl();

	void set_params(real* p_start, int p_id, real* p_cell_size);
};

#endif /* CELL_H_ */

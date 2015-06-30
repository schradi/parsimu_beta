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
	int timestep;			/**< current timestep*/

	Cell(real* p_cell_size, int p_id, real* p_start, int p_timestep);

	Cell(real* p_cell_size);

	Cell(){num_part=0;adding=NULL;pl=NULL;};

	/**
	 * Insert a new Particle in the Cell with specific datas
	 */
	void insertParticle(int p_id, real p_m, real* p_X, real* p_V);

	/**
	 * Copy a existing Particle in the cell
	 * @param p_p pointer to existin particle
	 */
	void insertParticle(Particle* p_p);

	void insertParticle(ParticleList* p_pl);

	/**
	 * Move all particles stored in adding to pl
	 * Adding is empty afterwards
	 */
	void adding_to_pl();

	/**
	 * convert a ParticleList in an array of double
	 * @param cod_pl the resulting array will be stored at the end of this vector
	 * @param pos current length of the vector
	 */
	void code_pl(real* cod_pl, int pos, int size);

	/**
	 * convert an array of double in an Particle
	 * @param p storage, where the Particle will be stored
	 * @param cod_pl converted array o double
	 * @param pos current starting-position of conversion
	 */
	void uncodePl(Particle* p, real* cod_pl, int pos, int size);

	/**
	 * deletes a particle out o pl with an specific id
	 * @param p_id id of the particle to be deleted
	 */
	void deleteParticle(int p_id); /**< deleting Particle with id=p_id (does nothing, if it doesn't contain it)*/

	/**
	 * deletes a ParticleList Item
	 * @param p_pre previous ParticleList-Item
	 * @param p_del ParticleList-Item to be deleted
	 */
	void deleteParticle(ParticleList* p_pre, ParticleList* p_del);

	/**
	 * deletes the whole ParicleList pl
	 */
	void deletePl();

	void set_params(real* p_start, int p_id, real p_cell_size);
};

#endif /* CELL_H_ */

/*
 * particle.h
 *
 *  Created on: 24.03.2015
 *      Author: jonas
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "defines.h"

struct Particle{
public:
	int  id;			/**< ID of the Particle to track the way*/
	real m;				/**< Mass of the particle*/
	real X[DIM];		/**< Position of the particle*/
	real V[DIM];		/**< Velocity of the particle*/
	real F[DIM];		/**< Force-Vector of the particle in the next timestep*/
	real F_old[DIM];	/**< Force-Vector of the particle in the current timestep*/
};

struct ParticleList{
	Particle* p;		/**< Paricle*/
	ParticleList* next; /**< next ParticleList Item*/
	ParticleList(){p=NULL; next=NULL;};
};



#endif /* PARTICLE_H_ */

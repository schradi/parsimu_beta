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
	long int  id;			/**< ID of the Particle to track the way*/
	real m;				/**< Mass of the particle*/
	real X[DIM];		/**< Position of the particle*/
	real V[DIM];		/**< Velocity of the particle*/
	real a[DIM];		/**< Force-Vector of the particle in the next timestep*/
	real a_old[DIM];
	real Epot;
	real Ekin;
	Particle(){
		id=-1;
		m=1;
		for(int d=0; d<DIM; d++){
			X[d]=0;
			V[d]=0;
			a[d]=0;
			a_old[d]=0;
			Ekin=0;
			Epot=0;
		}
	}
};

struct ParticleList{
	Particle* p;		/**< Paricle*/
	ParticleList* next; /**< next ParticleList Item*/
	ParticleList(){p=NULL; next=NULL;};
};



#endif /* PARTICLE_H_ */

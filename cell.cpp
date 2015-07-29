/*
 * cell.cpp
 *
 *  Created on: 24.03.2015
 *      Author: jonas
 */

#include "cell.h"

Cell::Cell(real* p_cell_size, int p_id, real* p_start){
	Cell();
	set_params(p_start, p_id, p_cell_size);
}

void Cell::adding_to_pl(){
	if(pl!=NULL){
		ParticleList* tmp_pl;
		tmp_pl=pl;
		while(tmp_pl->next!=NULL){
			tmp_pl=tmp_pl->next;
		}
		tmp_pl->next=adding;
	}else{
		pl=adding;
	}
	adding=NULL;
}

void Cell::insertParticle(Particle* p_p){
	ParticleList* tmp_pl = new ParticleList;
	tmp_pl->p=p_p;
	tmp_pl->next=pl;
	pl=tmp_pl;
	num_part++;
}

void Cell::insertParticle(ParticleList* p_pl){
	p_pl->next=pl;
	pl=p_pl;
	num_part++;
}

void Cell::deletePl(){
	ParticleList* tmp;
	while(pl!=NULL){
		tmp=pl;
		pl=pl->next;
		num_part--;
		delete tmp;
	}
}

void Cell::set_params(real* p_start, int p_id, real* p_cell_size){
	id=p_id;
	for(int d=0; d<DIM; d++){
		start[d] = p_start[d];
		cell_size[d] = p_cell_size[d];
	}
}







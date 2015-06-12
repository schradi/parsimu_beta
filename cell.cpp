/*
 * cell.cpp
 *
 *  Created on: 24.03.2015
 *      Author: jonas
 */

#include "cell.h"

void Cell::set_params(real* p_start, int p_id, real p_cell_size){
	id=p_id;
	for(int d=0; d<DIM; d++){
		start[d] = p_start[d];
		cell_size[d] = p_cell_size;
	}
	num_part=0;
	timestep=0;
	pl=NULL;
	adding=NULL;
}

void Cell::insertParticle(int p_id, real p_m, real* p_X, real* p_V){
	ParticleList* tmp_pl = new ParticleList;
	tmp_pl->p = new Particle;
	tmp_pl->p->m=p_m;
	tmp_pl->p->id=p_id;
	for(int d=0; d<DIM; d++){
		tmp_pl->p->X[d]=p_X[d];
		tmp_pl->p->V[d]=p_V[d];
	}
	tmp_pl->next=pl;
	pl=tmp_pl;
}

void Cell::insertParticle(Particle* p_p){
	ParticleList* tmp_pl = new ParticleList;
	tmp_pl->p=p_p;
	tmp_pl->next=pl;
	pl=tmp_pl;
}

void Cell::code_pl(real* cod_pl, int pos, int size){
	ParticleList* tmp;
	tmp=pl;
	while(tmp!=NULL){
		cod_pl[pos]=(real) tmp->p->id;
		if(size>1) cod_pl[pos+1]=tmp->p->X[0];
		if(size>2) cod_pl[pos+2]=tmp->p->X[1];
		if(size>3) cod_pl[pos+3]=tmp->p->V[0];
		if(size>4) cod_pl[pos+4]=tmp->p->V[1];
		if(size>5) cod_pl[pos+5]=tmp->p->m;
		pos+=size;
	}
}

void Cell::uncodePl(Particle* p, real* cod_pl, int pos, int size){
	p->id=(int) cod_pl[pos];
	if(size>1) p->X[0]=cod_pl[pos+1];
	if(size>2) p->X[1]=cod_pl[pos+2];
	if(size>3) p->V[0]=cod_pl[pos+3];
	if(size>4) p->V[1]=cod_pl[pos+4];
	if(size>5) p->m=cod_pl[pos+5];
}

void Cell::deletePl(){
	ParticleList* tmp;
	while(pl!=NULL){
		tmp=pl;
		pl=pl->next;
		delete tmp;
	}
}

void Cell::deleteParticle(int p_id){
	ParticleList* tmp_pl_pre = pl;
	ParticleList* tmp_pl_del = pl;
	while (tmp_pl_pre->next->p->id!=p_id){
		tmp_pl_pre=tmp_pl_pre->next;
	}
	tmp_pl_del=tmp_pl_pre->next;
	tmp_pl_pre->next=tmp_pl_pre->next->next;
	delete(tmp_pl_del);
}

void Cell::deleteParticle(ParticleList* p_pre, ParticleList* p_del){
	if(p_pre->next==p_del){
		p_pre->next=p_del->next;
		delete(p_del);
	}
}

void Cell::adding_to_pl(){
	ParticleList* tmp_pl;
	tmp_pl=pl;
	while(tmp_pl->next!=NULL){
		tmp_pl=tmp_pl->next;
	}
	tmp_pl->next=adding;
	adding=NULL;
}



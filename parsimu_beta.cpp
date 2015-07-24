//============================================================================
// Name        : parsimu_beta.cpp
// Author      : Jonas Schradi
// Version     :
// Copyright   : Jonas Schradi
// Description : Particle Simulation
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <mpi.h>

//#include <direct.h>

#include "defines.h"
#include "particle.h"
#include "cell.h"
#include "SimProcess.h"
#include "timers.h"

int main(int argc, char* argv[]) {



	int c_nump;
	if(argc>=2){
		c_nump=atoi(argv[1]);
	}else{
		c_nump=0;
	}
	MPI::Init (argc, argv);

	std::fstream del;
	del.open("data/energy.csv", std::ios::trunc);
	del.close();

	// calculate global_np
	int global_np[DIM];
	int pr=MPI::COMM_WORLD.Get_size();
	for(int d=0; d<DIM; d++){
		global_np[d]=pow(pr,(double) 1/DIM);
		if((pr/pow(global_np[d],DIM))!=1){
			std::cout<<"\n-----ABORTED!-----\nwrong number of processes.\n\n";
			MPI::Finalize();
			return 0;
		}
	}

	SimProcess* sim_p;
	sim_p= new SimProcess();
	sim_p->errt=0;
	int numc=(sim_p->local_nc[0]+2)*(sim_p->local_nc[1]+2);
	Cell cells[numc];
	if(sim_p->rank==0) std::cout<<"create_cells\n";
	sim_p->create_cells(cells);
	MPI::COMM_WORLD.Barrier();

	sim_p->num_part=c_nump;
	if(sim_p->rank==0) std::cout<<"initData\n";
	sim_p->initData(cells);

	sim_p->output(cells, 0);
	MPI::COMM_WORLD.Barrier();
	sim_p->communicate(cells);
	MPI::COMM_WORLD.Barrier();
	sim_p->timeIntegration(cells);
//	MPI::COMM_WORLD.Barrier();
//	if(sim_p->rank==0){
//		std::cout<<"DONE\n";
//		std::fstream file;
//		file.open("times.csv", std::ios::out | std::ios::app);
//		file<<c_nump*c_nump;
//		for(TimerList* ti=sim_p->timerList; ti!=NULL; ti=ti->next){
//			file<<" "<<ti->t->timer;
//			std::cout<<ti->t->tag<<" ";
//		}
//		file<<"\n";
//		std::cout<<"\n";
//		file.close();
//	}
//
	MPI::Finalize();
	return 0;
}

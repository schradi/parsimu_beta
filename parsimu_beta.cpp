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

#include "defines.h"
#include "particle.h"
#include "cell.h"
#include "SimProcess.h"

int main(int argc, char* argv[]) {
	MPI::Init (argc, argv);

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
		if(global_np[d]%2!=0){
			//Schachbrettmuster funktioniert sonst nicht, wird bei kommunikation benötigt
			std::cout<<"\n-----ABORTED!-----\nwrong number of processes.\n\n";
			MPI::Finalize();
			return 0;
		}
	}
	real global_size[DIM];
	for(int d=0; d<DIM; d++){
		global_size[d]=100;
	}
	real r_cut=25;
	SimProcess* sim_p;

	sim_p= new SimProcess(r_cut, global_size, global_np);
	int numc=(sim_p->local_nc[0]+2)*(sim_p->local_nc[1]+2);
	Cell cells[numc];
	sim_p->create_cells(cells);
	MPI::COMM_WORLD.Barrier();
	sim_p->initData(cells);
	real delta_t=1;
	real t_end=1000;
	sim_p->t=0;
	sim_p->output_resolution=5;
	sim_p->output(delta_t, cells);
	MPI::COMM_WORLD.Barrier();
	sim_p->communicate(cells);

	MPI::COMM_WORLD.Barrier();

	sim_p->timeIntegration(delta_t, t_end, cells);

	MPI::Finalize();
	return 0;
}

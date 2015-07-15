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
		if(global_np[d]%2!=0 && global_np[d]!=1){
			//Schachbrettmuster funktioniert sonst nicht, wird bei kommunikation benötigt
			std::cout<<"\n-----ABORTED!-----\nwrong number of processes.\n\n";
			MPI::Finalize();
			return 0;
		}
	}

	SimProcess* sim_p;
	sim_p= new SimProcess();
	int numc=(sim_p->local_nc[0]+2)*(sim_p->local_nc[1]+2);
	Cell cells[numc];
	sim_p->create_cells(cells);
	MPI::COMM_WORLD.Barrier();
	sim_p->num_part=c_nump;
	sim_p->initData(cells);

	sim_p->output(cells, 0);
	MPI::COMM_WORLD.Barrier();
	if(global_np[0]!=1) sim_p->communicate(cells);
	sim_p->output(cells, 1);

	sim_p->timeIntegration(cells);
	MPI::COMM_WORLD.Barrier();
	if(sim_p->rank==0){
		std::cout<<"times.csv\n";
		std::fstream file;
//			if(npart==sqrt(c_start)){
//				file.open("times.csv", std::ios::out | std::ios::trunc);
//				file<<"Number Particles";
//				for(TimerList* ti=sim_p->timerList; ti!=NULL; ti=ti->next){
//					file<<","<<ti->t->tag;
//				}
//				file<<"\n";
//			}else{
			file.open("times.csv", std::ios::out | std::ios::app);
//			}
			file<<c_nump*c_nump;
		for(TimerList* ti=sim_p->timerList; ti!=NULL; ti=ti->next){
			file<<" "<<ti->t->timer;
			std::cout<<ti->t->tag<<" ";
		}
		file<<"\n";
		std::cout<<"\n";
		file.close();
	}

	MPI::Finalize();
//	std::cout<<"C_Start\t"<<c_start<<"\n";
//	std::cout<<"C_Res\t"<<c_res<<"\n";
//	std::cout<<"C_Stop\t"<<c_stop<<"\n";
	return 0;
}

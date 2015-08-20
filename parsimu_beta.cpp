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
	long int c_nump=0;;
	char* c_folder=new char[20];
	bool c_log_time=false;
	bool c_log_energy=false;
	bool c_log_positions=false;
	bool c_log_velocity=false;
	bool c_log_id=false;
	int i=1;
	while (i<argc-1){
		if(strcmp(argv[i], "-p")==0){
			c_nump=atoi(argv[i+1]);
		}
		if(strcmp(argv[i], "-output")==0){
			strcpy(c_folder,argv[i+1]);
			c_log_positions=true;
		}
		if(strcmp(argv[i], "-log")==0){
			if(strcmp(argv[i+1], "time")==0){
				c_log_time=true;
			}
			if(strcmp(argv[i+1], "energy")==0){
				c_log_energy=true;
				if(c_folder==NULL){
					strcpy(c_folder, "data");
				}
			}
			if(strcmp(argv[i+1], "v")==0){
				c_log_velocity=true;
			}
			if(strcmp(argv[i+1], "id")==0){
				c_log_id=true;
			}
		}
		i+=2;
	}
	MPI::Init (argc, argv);
	if(c_log_energy){
		std::fstream del;
		char filename[20];
		strcpy(filename, c_folder);
		strcat(filename, "energy.csv");
		del.open(filename, std::ios::out | std::ios::trunc);
		del.close();
	}

	SimProcess* sim_p;
	sim_p= new SimProcess(c_folder);
	// setting parameters
	sim_p->FALSCH=false;
	sim_p->log_time=c_log_time;
	sim_p->log_energy=c_log_energy;
	sim_p->log_positions=c_log_positions;
	sim_p->log_velocity=c_log_velocity;
	sim_p->log_id=c_log_id;
	if(sim_p->rank==0)	sim_p->num_part=c_nump;
	sim_p->global_num_part=pow(c_nump, DIM);


	real* p_map=NULL;
	if(sim_p->rank==0){
		p_map=new real [DIM*sim_p->np];
	}
	MPI::COMM_WORLD.Barrier();
	sim_p->devide_symetric(p_map);
	sim_p->spread_local_info(p_map);
	if(sim_p->rank==0) std::cout<<"create_cells\n";
	int numc=(sim_p->local_nc[0]+2)*(sim_p->local_nc[1]+2);
	Cell* cells=new Cell[numc];
	sim_p->create_cells(cells);
	MPI::COMM_WORLD.Barrier();

	if(sim_p->rank==0) std::cout<<"initData\n";
	sim_p->initData(cells, p_map);

	sim_p->output(cells, 0);
	MPI::COMM_WORLD.Barrier();
	if(sim_p->rank==0) std::cout<<"communicate\n";

	sim_p->communicate(cells);

	sim_p->timeIntegration(cells, p_map);
//	int c=0;
//	std::cout<<	std::cout<<"P"<<rank<<" delta_t "<<delta_t<<"\n";"done";
//	for(long int i=0; i<200000000; i+=1){
//		i+=c;
//	}
//	MPI::COMM_WORLD.Barrier();
//	if(sim_p->rank==0){
//		std::cout<<"DONE\n";
//		std::fstream file;
//		file.open("Times.csv", std::ios::out | std::ios::app);
//		file<<c_nump*c_nump;
//		for(TimerList* ti=sim_p->timer_list; ti!=NULL; ti=ti->next){
//			file<<ti->t->timer<<",\t";
//			std::cout<<ti->t->tag<<", ";
//		}
//		file<<"\n";
//		std::cout<<"\n";
//		file.close();
//	}

	MPI::Finalize();
	return 0;
}

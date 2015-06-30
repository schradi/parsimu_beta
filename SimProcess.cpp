/*
 * SImProcess.cpp
 *
 *  Created on: 21.04.2015
 *      Author: jonas
 */

#include "SimProcess.h"

void SimProcess::timeIntegration(real delta_t, real t_end, Cell* cells){
	// todo: Zeit stoppen, die zur Berechnung benötigt wird, zur Lastbalancierung interessant

	t=0;
	max_V=cell_size[0]/delta_t;
	while (t<t_end){
		MPI::COMM_WORLD.Barrier();
//		compF(cells);
//		compV(delta_t, cells);
//		std::cout<<"Pr "<<rank<<" - compX\n";
		compX(delta_t, cells);
//		std::cout<<"Pr "<<rank<<" - moveParticles\n";
		moveParticles(cells);
		std::cout<<"Pr "<<rank<<" - communicate\n";
		communicate(cells);
		if(((int) (t/delta_t))%output_resolution==0){
			std::cout<<"Pr "<<rank<<" - output\n";
			output(delta_t, cells);
		}
		t+=delta_t;
//		status_output("processing");
	}
}

void SimProcess::output(real delta_t, Cell* cells){
	MPI::COMM_WORLD.Barrier();
	int ic[DIM];

//	// Important: just particles stored in pl are considered

	MPI::Request request;
	MPI::Status status;

	int* send_pl_l = new int;
	*send_pl_l=get_num_p(ic_start, ic_stop, cells);
		/* direction of communication is not important,
		 * just needs to be different from other direction
		 */
//	std::cout<<"Process "<<rank<<" contains "<<*send_pl_l<<" particles\n";
	std::cout<<"Pr "<<rank<<" - output - begin\n";
	if(rank==0){

		// MASTER: Receives information from other processes and puts them in a file
		int recv_pl_l[global_np[0]*global_np[1]]; // Number of particles of each process
		int recv_pl_al=*send_pl_l;	// Number of particles at all
		recv_pl_l[0]=*send_pl_l;

		// Receiving Number of Particles of each Process
		int* new_recv;
		new_recv=new int;
		for (int cp=1; cp<MPI::COMM_WORLD.Get_size();cp++){
			*new_recv=0;
			request=MPI::COMM_WORLD.Irecv(new_recv, 1, MPI::INT, cp, 1);
			request.Wait(status);
			recv_pl_l[cp]=*new_recv;
			recv_pl_al+=*new_recv;
		}
		delete new_recv;
		// Receiving Particlelists
		real recv_pl[OUTP_SZE*recv_pl_al];
		int ic[DIM];
		int pos=0;
		for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
			for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
				for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
					if (OUTP_SZE>0) recv_pl[pos]=i->p->id;
					if (OUTP_SZE>1) recv_pl[pos+1]=i->p->X[0];
					if (OUTP_SZE>2) recv_pl[pos+2]=i->p->X[1];
					if (OUTP_SZE>3) recv_pl[pos+3]=i->p->V[0];
					if (OUTP_SZE>4) recv_pl[pos+4]=i->p->V[1];
					pos+=OUTP_SZE;
				}
			}
		}

		for (int cp=1; cp<global_np[0]*global_np[1];cp++){
			request=MPI::COMM_WORLD.Irecv(&recv_pl[pos], recv_pl_l[cp]*OUTP_SZE, MPI::DOUBLE, cp, 18);
			request.Wait(status);
			pos+=recv_pl_l[cp]*OUTP_SZE;
		}


		// Fileoutput
		std::fstream file;
		int outp_nr=t/(delta_t*output_resolution);

		// Geting Name of the file
		char cline_nr[32];
		sprintf(cline_nr, "%d", outp_nr+1);
		char outputfile_new[20]="data/data.csv.";
		strcat(outputfile_new, cline_nr);

		file.open(outputfile_new, std::ios::out|std::ios::trunc);
		file<<"x coord,y coord\n";
		for(pos=0; pos<recv_pl_al*OUTP_SZE; pos+=OUTP_SZE){
//			cells[0].uncodePl(p, recv_pl, pos, OUTP_SZE);
			file<<recv_pl[pos+1];
			for(int i=2; i<OUTP_SZE; i++){
				file<<","<<recv_pl[pos+i];
			}
			file<<"\n";
		}
		file.close();

	}else{
		// Sending number of Particles
		real send_pl[*send_pl_l*OUTP_SZE];
		std::cout<<"Pr "<<rank<<" - output - sending: "<<*send_pl_l<<" Particles\n";
		request=MPI::COMM_WORLD.Isend(send_pl_l, 1, MPI::INT, 0, 1);
		request.Wait(status);

		// Sending Particles
		int ic[DIM];
		int pos=0;
		for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
			for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
				for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
//					std::cout<<"Pr "<<rank<<" - output - sending - inner loop - begin\n";
					if (OUTP_SZE>0) send_pl[pos]=i->p->id;
					if (OUTP_SZE>1) send_pl[pos+1]=i->p->X[0];
					if (OUTP_SZE>2) send_pl[pos+2]=i->p->X[1];
					if (OUTP_SZE>3) send_pl[pos+3]=i->p->V[0];
					if (OUTP_SZE>4) send_pl[pos+4]=i->p->V[1];
					pos+=OUTP_SZE;
//					std::cout<<"Pr "<<rank<<" - output - sending - outer loop - end\n";
				}
			}
		}
		request=MPI::COMM_WORLD.Isend(send_pl, *send_pl_l*OUTP_SZE, MPI::DOUBLE, 0, 18);
		request.Wait(status);
	}
	delete send_pl_l;
	std::cout<<"Pr "<<rank<<" - output - end\n";
}

void SimProcess::compF(Cell* cells){
	int ic[DIM];
	int nc[DIM];
	Cell* ci;
	Cell* cj;
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			ci=&cells[local_index(ic)];
			for(ParticleList* pi=ci->pl; pi->next!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for(int d=0; d<DIM; d++){
					pi->p->F_old[d]=pi->p->F[d];
					pi->p->F[d]=0;
				}
			}
		}
	}
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			ci=&cells[local_index(ic)];

			for(ParticleList* pi=ci->pl; pi->next!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for (nc[1]=ic[1]-1;nc[1]<=ic[1]+1;nc[1]++){
					for (nc[0]=ic[0]-1;nc[0]<=ic[0]+1;nc[0]++){
						// for each neighboring cell
						cj=&cells[local_index(nc)];
						for(ParticleList* pj=cj->pl; pj!=NULL; pj=pj->next ){
							compF(pi->p, pj->p);
						}
					}
				}
			}
		}
	}
}

void SimProcess::compF(Particle* p_i, Particle* p_j){
	real F[DIM];
	F[0]=p_i->X[0]-p_j->X[0];
	F[1]=p_i->X[1]-p_j->X[1];
	force(F);
	for(int d=0; d<DIM; d++){
		p_i->F[d]+=F[d];
		p_j->F[d]+=F[d];
	}
}

void SimProcess::compV(real delta_t, Cell* cells){
	Cell* c;
	int ic[DIM];
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
			for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			c=&cells[local_index(ic)];
			for(ParticleList* pi=c->pl; pi->next!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for(int d=0; d<DIM; d++){
					pi->p->V[d]=pi->p->V[d]+pi->p->F[d]*delta_t;
				}
			}
		}
	}
}

void SimProcess::compX(real delta_t, Cell* cells){
	Cell* c;
	int ic[DIM];
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			// for each Cell
			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for(int d=0; d<DIM; d++)
				pi->p->X[d]=pi->p->X[d]+pi->p->V[d]*delta_t;
			}
		}
	}
}

void SimProcess::moveParticles(Cell* cells){
	Cell* c;
	int   ic[DIM];
	int   new_ic[DIM];
	bool  new_cell;
	ParticleList* tmp;
	ParticleList* prev;
	ParticleList* akt;
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			prev=NULL;
			// for each Cell
			c=&cells[local_index(ic)];
			akt=c->pl;
			while(akt!=NULL){
				new_cell=false;
				for(int d=0; d<DIM; d++){
					new_ic[d]=ic[d];
					if(akt->p->X[d]<c->start[d]){
						new_ic[d]-=1;
						new_cell=true;
					}else{
						if(akt->p->X[d]>=c->start[d]+cell_size[d]){
							new_ic[d]+=1;
							new_cell=true;
						}
					}
				}
				if(new_cell){
					if(prev!=NULL){
						prev->next=akt->next;
					}else{
						c->pl=akt->next;
					}
					tmp=cells[local_index(new_ic)].adding;
					cells[local_index(new_ic)].adding=akt;
					cells[local_index(new_ic)].adding->next=tmp;
					cells[local_index(new_ic)].num_part++;
					c->num_part--;
					if(prev!=NULL){
						akt=prev->next;
					}else{
						akt=c->pl;
					}
				}else{
					prev=akt;
					akt=akt->next;
				}
			}
		}
	}
}

int SimProcess::getLocalCellId(real* X){
	int id[DIM];
	id[0]=(int) ceil((X[0]/r_cut)-1);	// Abrunden
	id[1]=(int) ceil((X[1]/r_cut)-1);
	return local_index(id);
}

SimProcess::SimProcess(real p_r_cut, real* p_global_size, int* p_global_np){
	rank = MPI::COMM_WORLD.Get_rank();
	num_part=0;
	r_cut=p_r_cut;

	// PositionsUNabhängige Parameter berechnen
	for(int d=0; d<DIM; d++){
		global_size[d]=p_global_size[d];
		global_np[d]=p_global_np[d];
		local_size[d]=global_size[d]/global_np[d];

		cell_size[d]=r_cut;
		real pos=0;
		local_nc[d]=0;
		while(pos+cell_size[d]<=local_size[d]){
			local_nc[d]++;
			pos+=cell_size[d];
		}
		cell_size[d]=local_size[d]/local_nc[d];
		global_nc[d]=local_nc[d]*global_np[d];
		ip[d]=0;
	}
	np=global_np[0]*global_np[1];
	// PositionsABhängige Parameter berechnen
	while(ip[0]+global_np[0]*ip[1]!=rank){
		if(ip[0]+1<global_np[0]){
			ip[0]++;
		}else{
			ip[0]=0;
			ip[1]++;
		}
	}

	for(int d=0; d<DIM; d++){
		start[d]=ip[d]*local_size[d];
		ic_start[d]=ip[d]*local_nc[d];
		ic_stop[d]=ic_start[d]+local_nc[d]-1;
	}

	// Nachbarn berechnen
	if(ip[0]-1>=0){
		neigh_lower[0]=ip[0]-1+global_np[0]*ip[1];
	}else{
		neigh_lower[0]=(global_np[0]*(ip[1]+1))-1;
	}
	if(ip[0]+1<global_np[0]){
		neigh_upper[0]=ip[0]+1+global_np[0]*ip[1];
	}else{
		neigh_upper[0]=(global_np[0]*ip[1]);
	}

	if(ip[1]-1>=0){
		neigh_lower[1]=ip[0]+global_np[0]*(ip[1]-1);
	}else{
		neigh_lower[1]=global_np[0]*(global_np[1]-1)+ip[0];
	}
	if(ip[1]+1<global_np[1]){
		neigh_upper[1]=ip[0]+global_np[0]*(ip[1]+1);
	}else{
		neigh_upper[1]=ip[0];
	}

	// Ausgabe zur Überprüfung möglich:
	if(rank==0){
		std::cout<<"-------RANK: "<<rank<<"-------\n";
		for(int d=0; d<DIM; d++){
			std::cout<<"neigh_lower["<<d<<"]\t"<<neigh_lower[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"neigh_upper["<<d<<"]\t"<<neigh_upper[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"start["<<d<<"]\t"<<start[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"ic_start["<<d<<"]\t"<<ic_start[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"ic_stop["<<d<<"]\t"<<ic_stop[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"cell_size["<<d<<"]\t"<<cell_size[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"global_nc["<<d<<"]\t"<<global_nc[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"local_nc["<<d<<"]\t"<<local_nc[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"local_size["<<d<<"]\t"<<local_size[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			std::cout<<"ip["<<d<<"]\t\t"<<ip[d]<<"\n";
		}
		std::cout<<"id of last cell:\t"<<index(ic_stop, global_nc)<<"\n";
	}
//	status_output((char*)"ready to start");
}

void SimProcess::create_cells(Cell* cells){
	// Create Cells
	real p_cell_start[DIM];
	int ic[DIM];
	for (ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for (ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			for (int d=0; d<DIM; d++){
				p_cell_start[d]=ic[d]*r_cut;
			}
//			std::cout<<"index("<<ic[0]<<"/"<<ic[1]<<")="<<index(ic, global_nc)<<"\n";
			cells[local_index(ic)].set_params(p_cell_start, index(ic, global_nc), r_cut);
		}
	}
}

void SimProcess::force(real* X){
	real sigma=1;	/**< Parameter sigma for the Lennard-Jones-Potential*/
	real epsilon=5; /**< Parameter epsilon for the Lennard-Jones-Potential*/
	real r=0; /**< length*/
	for (int d=0;d<DIM;d++){
		r+=X[d]*X[d];
	}
	r=sqrt(r);

	// Normalization of X
	for (int d=0; d<DIM; d++){
		X[d]=X[d]/r;
	}

	real s = sqr(sigma) /r;
	s=sqr(s)*s;
	real f = 24* epsilon * s/ r * (1-2*s);
	for (int d=0; d<DIM; d++){
		X[d] = f*X[d];
	}
}

void SimProcess::communicate(Cell* cells){
	int ic[DIM];
	// For all surrounding Ghost-Cells: delete Pl
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			if(ic[1]==ic_start[1]-1||ic[1]==ic_stop[1]+1
					||ic[0]==ic_start[0]-1||ic[0]==ic_stop[0]+1){
//				if(rank==0) std::cout<<"Index: "<<local_index(ic)<<"\n";
				cells[local_index(ic)].deletePl();
			}
		}
	}

	int icr_start[DIM]={ic_start[0]-1, ic_start[1]-1};
	int icr_stop[DIM]={ic_stop[0]+1, ic_stop[1]+1};
	num_part=get_num_p(icr_start, icr_stop, cells);

	// For all cells: include adding in PL
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			adding_to_pl(cells);
		}
	}

	MPI::COMM_WORLD.Barrier();
	communicate(0, cells);
	MPI::COMM_WORLD.Barrier();
	communicate(1, cells);
}

void SimProcess::communicate(int com_d, Cell* cells){
	int oth_d;
	if(com_d==0){
		oth_d=1;
	}else if(com_d==1){
		oth_d=0;
	}
	int ic[DIM];
	int icr_start[DIM];
	int icr_stop[DIM];

	/**
	 * Schachbrettkommunikation nötig, da Prozesse nach senden aufeinander warten müssen, bis alles gesandt wurde.
	 */

	if(ip[com_d]%2==0){
		//right
		icr_start[com_d]=ic_stop[com_d];
		icr_start[oth_d]=ic_start[oth_d]-1;
		icr_stop[com_d]=ic_stop[com_d]+1;
		icr_stop[oth_d]=ic_stop[oth_d]+1;

		int send_plr_l = get_num_p(icr_start, icr_stop, cells);
		real send_plr[send_plr_l*CD_P_SZE];
		code_pl(send_plr, icr_start, icr_stop, CD_P_SZE, cells);
		pl_send_l(&send_plr_l, neigh_upper[com_d], 1);
		pl_send(send_plr, neigh_upper[com_d], send_plr_l, 2);

		delete_pl(icr_start, icr_stop, cells);

		int recv_plr_l;
		pl_recv_l(&recv_plr_l, neigh_upper[com_d], 3);
		real recv_plr[recv_plr_l*CD_P_SZE];
		pl_recv(recv_plr, recv_plr_l, neigh_upper[com_d], 4);
		uncodePl(recv_plr, icr_start, icr_stop, recv_plr_l, CD_P_SZE, cells);

		//left
		icr_start[com_d]=ic_start[com_d]-1;
		icr_start[oth_d]=ic_start[oth_d]-1;
		icr_stop[com_d]=ic_start[com_d];
		icr_stop[oth_d]=ic_stop[oth_d]+1;

		int send_pll_l = get_num_p(icr_start, icr_stop, cells);
		real send_pll[send_pll_l*CD_P_SZE];
		code_pl(send_pll, icr_start, icr_stop, CD_P_SZE, cells);
		pl_send_l(&send_pll_l, neigh_lower[com_d], 5);
		pl_send(send_pll, neigh_lower[com_d], send_pll_l, 6);

		delete_pl(icr_start, icr_stop, cells);

		int recv_pll_l;
		pl_recv_l(&recv_pll_l, neigh_lower[com_d], 7);
		real recv_pll[recv_pll_l*CD_P_SZE];
		pl_recv(recv_pll, recv_pll_l, neigh_lower[com_d], 8);
		uncodePl(recv_pll, icr_start, icr_stop, recv_pll_l, CD_P_SZE, cells);

	}else{
		// left
		icr_start[com_d]=ic_start[com_d]-1;
		icr_start[oth_d]=ic_start[oth_d]-1;
		icr_stop[com_d]=ic_start[com_d];
		icr_stop[oth_d]=ic_stop[oth_d]+1;

		int recv_pll_l;
		pl_recv_l(&recv_pll_l, neigh_lower[com_d], 1);
		real recv_pll[recv_pll_l*CD_P_SZE];
		pl_recv(recv_pll, recv_pll_l, neigh_lower[com_d], 2);
		uncodePl(recv_pll, icr_start, icr_stop, recv_pll_l, CD_P_SZE, cells);
		int send_pll_l = get_num_p(icr_start, icr_stop, cells);
		real send_pll[send_pll_l*CD_P_SZE];
		code_pl(send_pll, icr_start, icr_stop, CD_P_SZE, cells);
		pl_send_l(&send_pll_l, neigh_lower[com_d], 3);
		pl_send(send_pll, neigh_lower[com_d], send_pll_l, 4);


		// right
		icr_start[com_d]=ic_stop[com_d];
		icr_start[oth_d]=ic_start[oth_d]-1;
		icr_stop[com_d]=ic_stop[com_d]+1;
		icr_stop[oth_d]=ic_stop[oth_d]+1;

		int recv_plr_l;
		pl_recv_l(&recv_plr_l, neigh_lower[com_d], 5);
		real recv_plr[recv_plr_l*CD_P_SZE];
		pl_recv(recv_plr, recv_plr_l, neigh_lower[com_d], 6);
		uncodePl(recv_plr, icr_start, icr_stop, recv_plr_l, CD_P_SZE, cells);
		int send_plr_l = get_num_p(icr_start, icr_stop, cells);
		real send_plr[send_plr_l*CD_P_SZE];
		code_pl(send_plr, icr_start, icr_stop, CD_P_SZE, cells);
		pl_send_l(&send_plr_l, neigh_lower[com_d], 7);
		pl_send(send_plr, neigh_lower[com_d], send_plr_l, 8);
	}

	adding_to_pl(cells);
	// Definition: 	right = ascending in com_d direction
	// 				left = descending in com_d direction
}

void SimProcess::delete_pl(int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			cells[local_index(ic)].deletePl();
		}
	}
}

int SimProcess::get_num_p(int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	int num=0;
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			num+=cells[local_index(ic)].num_part;
		}
	}
	return num;
}

void SimProcess::pl_send_l(int* length, int recv_rank, int tag){
//	std::cout<<"Process "<<rank<<" is sending length "<<*length<<" to "<<recv_rank<<" with tag "<<tag<<"\n";
	MPI::Request request;
	MPI::Status status;
	request=MPI::COMM_WORLD.Isend(length, 1, MPI::INT, recv_rank, tag);
	request.Wait(status);
}

void SimProcess::pl_send(real* to_send, int recv_rank, int length, int tag){
//	std::cout<<"Process "<<rank<<" is sending list containing "<<length<<" particles to "<<recv_rank<<" with tag "<<tag<<"\n";
	MPI::Request request;
	MPI::Status status;
	request=MPI::COMM_WORLD.Isend(to_send, length*CD_P_SZE, MPI::DOUBLE, recv_rank, tag);
	request.Wait(status);
}

void SimProcess::pl_recv_l(int* length, int send_rank, int tag){
//	std::cout<<"Process "<<rank<<" tries to receive length from "<<send_rank<<" with tag "<<tag<<"\n";
	MPI::Request request;
	MPI::Status status;
	request=MPI::COMM_WORLD.Irecv(length, 1, MPI::INT, send_rank, tag);
	request.Wait(status);
}

void SimProcess::pl_recv(real* recv, int length, int send_rank, int tag){
	MPI::Request request;
	MPI::Status status;
	request=MPI::COMM_WORLD.Irecv(recv, length*CD_P_SZE, MPI::DOUBLE, send_rank, tag);
	request.Wait(status);
//	std::cout<<"Process "<<rank<<" received list containing "<<length<<" particles from "<<send_rank<<" with tag "<<tag<<"\n";
}

//void SimProcess::pl_l_send_recv(int* send_pl_l, int* recv_pl_l, int com_d){
//	MPI::Request request;
//	MPI::Status status;
//
//	// send length of the collected List
//	request=MPI::COMM_WORLD.Isend(send_pl_l, 1, MPI::INT, neigh_upper[com_d], 1);
//
//	// receive length of the received list
//	request=MPI::COMM_WORLD.Irecv(recv_pl_l, 1, MPI::INT, neigh_lower[com_d], 1);
//	request.Wait(status);
//}

//void SimProcess::pl_send_recv(real* send_pl, int send_pl_l, real* recv_pl, int recv_pl_l, int com_d){
//	MPI::Request request;
//	MPI::Status status;
//	// send pl coded as array of double to the right
//	request=MPI::COMM_WORLD.Isend(send_pl, send_pl_l, MPI::DOUBLE, neigh_upper[com_d], 2);
//
//	// receive pl coded as array of double from the leftrsion from ‘int’ to ‘int*’ [-fpermissive]
//
//	request=MPI::COMM_WORLD.Irecv(recv_pl, recv_pl_l, MPI::DOUBLE, neigh_lower[com_d], 2);
//
//	// wait till communication is completed
//	request.Wait(status);
//}

void SimProcess::uncodePl(real* recv_pl, int* icr_start, int*icr_stop, int length_recv, int size, Cell* cells){
	int pos=0;
	Particle* p;
	int ic[DIM];
	while(pos<length_recv*CD_P_SZE){
		p=new Particle();
		num_part++;
		if(size>0) p->id=recv_pl[pos];
		if(size>1) p->X[0]=recv_pl[pos+1];
		if(size>2) p->X[1]=recv_pl[pos+2];
		if(size>3) p->V[0]=recv_pl[pos+3];
		if(size>4) p->V[1]=recv_pl[pos+4];
		if(size>5) p->m=recv_pl[pos+5];


		//Periodic Boundaries
		for(int d=0; d<DIM; d++){
			if(ic_start[d]==0&&p->X[d]>global_size[d]-cell_size[d]){
				p->X[d]-=global_size[d];
			}
			if(neigh_upper[d]<rank && p->X[d]<cell_size[d]){
				p->X[d]+=global_size[d];
			}
		}

		ic[0]=0;
		ic[1]=0;
		while(p->X[0]>=cell_size[0]*ic[0]){
			ic[0]++;
		}
		ic[0]-=1;
		while(p->X[1]>=cell_size[1]*ic[1]){
			ic[1]++;
		}
		ic[1]-=1;

		cells[local_index(ic)].insertParticle(p);
//		std::cout<<"Process "<<rank<<" included a Particle with Position ("<<p->X[0]<<"/"<<p->X[1]<<") in cell "<<local_index(ic)<<"\n";
		pos+=size;
	}
}

void SimProcess::code_pl(real* send_pl, int* icr_start, int* icr_stop, int size, Cell* cells){
	int ic[DIM];
	int pos=0;
	Cell* c;
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
				if(size>0) send_pl[pos]=pi->p->id;
				if(size>1) send_pl[pos+1]=pi->p->X[0];
				if(size>2) send_pl[pos+2]=pi->p->X[1];
				if(size>3) send_pl[pos+3]=pi->p->V[0];
				if(size>4) send_pl[pos+4]=pi->p->V[1];
				if(size>5) send_pl[pos+5]=pi->p->m;
				pos+=size;
			}
		}
	}
}

void SimProcess::create_particles(ParticleList* new_pl, real* r_start, real* r_stop, real resolution, real* p_V){
	real pos[DIM];
//	new_pl->next=NULL;
	if(new_pl->next==NULL) std::cout<<"jep\n";
	ParticleList* tmp;
	int id_c=0;
	for(pos[1]=r_start[1]; pos[1]<r_stop[1] && pos[1]<global_size[1]; pos[1]+=resolution){
		for(pos[0]=r_start[0]; pos[0]<r_stop[0] && pos[0]<global_size[0]; pos[0]+=resolution){
			// Da new_pl auf eine feste Adresse zeigt, kann diese nicht verändert werden. Lösung: Arbeiten mit new_pl->next
			if(new_pl->p!=NULL){
				tmp=new_pl->next;
				new_pl->next=new ParticleList;
				new_pl->next->next=tmp;
				new_pl->next->p=new_pl->p;
				new_pl->p=NULL;
			}
			new_pl->p= new Particle;
			new_pl->p->id=id_c;
			id_c++;
			new_pl->p->m=1;
			for(int d=0; d<DIM; d++){
				new_pl->p->X[d]=pos[d];
				new_pl->p->V[d]=p_V[d];
			}
		}
	}
	std::cout<<"Included Particles\t"<<id_c<<"\n";
}

void SimProcess::initData(Cell* cells){
	MPI::Request request;
	MPI::Status status;
	ParticleList* own=NULL;

	if(rank==0){
		// MASTER creates ParticleList containing all needed information
		real r_start[DIM]={1.5, 51.5};
		real r_stop[DIM]={24.9, 94.9};
		real resolution=1;
		real V[DIM]={1, 0.5};
		ParticleList* to_insert=new ParticleList;
		create_particles(to_insert, r_start, r_stop, resolution, V);

		// Sort Particles in different ParticleLists
		ParticleList* tmp;
		ParticleList* sorted_Particles[global_np[0]*global_np[1]];
		int *pl_l[global_np[0]*global_np[1]];
		for(int p_rank=0; p_rank<global_np[0]*global_np[1]; p_rank++){
			pl_l[p_rank]=new int;
			*pl_l[p_rank]=0;
			sorted_Particles[p_rank]=NULL;
		}
		int proc[DIM];
		int proc_idx=0;
		int n=0;
		while(to_insert!=NULL){
			n++;
			proc[0]=0;
			proc[1]=0;
			while(to_insert->p->X[0]>=(proc[0]+1)*local_size[0]){
				proc[0]++;
			}
			while(to_insert->p->X[1]>=(proc[1]+1)*local_size[1]){
				proc[1]++;
			}
			proc_idx=index(proc, global_np);
			// current particle is placed in the beginning of the list for each process
			tmp=to_insert->next;
			to_insert->next=sorted_Particles[proc_idx];
			sorted_Particles[proc_idx]=to_insert;
			to_insert=tmp;
			*pl_l[proc_idx]+=1;
		}
		// Particles are send to the different Processes
		// Sending number of Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			request=MPI::COMM_WORLD.Isend(pl_l[p], 1, MPI::INT, p, p);
			request.Wait(status);
		}
		// Sending Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			real cd_pl[*pl_l[p]*CD_P_SZE];
			ParticleList* tmp;
			tmp=sorted_Particles[p];
			for(int i=0; i<*pl_l[p]*CD_P_SZE; i+=CD_P_SZE){
				cd_pl[i]=(real) tmp->p->id;
				if(CD_P_SZE>1) cd_pl[1+i]=tmp->p->X[0];
				if(CD_P_SZE>2) cd_pl[2+i]=tmp->p->X[1];
				if(CD_P_SZE>3) cd_pl[3+i]=tmp->p->V[0];
				if(CD_P_SZE>4) cd_pl[4+i]=tmp->p->V[1];
				if(CD_P_SZE>5) cd_pl[5+i]=tmp->p->m;
				tmp=tmp->next;
			}
			request=MPI::COMM_WORLD.Isend(cd_pl, *pl_l[p]*CD_P_SZE, MPI::DOUBLE, p, global_nc[0]*global_nc[1]+p);
			request.Wait(status);
		}
		num_part=*pl_l[0];
		own=sorted_Particles[0];
	}else{
		// Receiving number of Particles
		int* new_recv;
		new_recv=new int;
		*new_recv=0;
		request=MPI::COMM_WORLD.Irecv(new_recv, 1, MPI::INT, 0, rank);
		request.Wait(status);
		num_part=*new_recv;
		delete new_recv;

		// Receiving Particles as Array of real
		real cod_pl[CD_P_SZE*num_part];
		request=MPI::COMM_WORLD.Irecv(cod_pl, CD_P_SZE*num_part, MPI::DOUBLE, 0, global_nc[0]*global_nc[1]+rank);
		request.Wait(status);

		// decode to ParticleList
		int ic[DIM];
		ParticleList* tmp;
		for(int pos=0; pos<CD_P_SZE*num_part; pos+=CD_P_SZE){
			tmp=own;
			own=new ParticleList;
			own->next=tmp;
			own->p=new Particle;
			if(CD_P_SZE>0) own->p->id=cod_pl[pos];
			if(CD_P_SZE>1) own->p->X[0]=cod_pl[pos+1];
			if(CD_P_SZE>2) own->p->X[1]=cod_pl[pos+2];
			if(CD_P_SZE>3) own->p->V[0]=cod_pl[pos+3];
			if(CD_P_SZE>4) own->p->V[1]=cod_pl[pos+4];
			if(CD_P_SZE>5) own->p->m=cod_pl[pos+5];
		}
	}

	// sort particles in Cells
	ParticleList* tmp;
	int ic[DIM];
	while(own!=NULL){
		ic[0]=0; ic[1]=0;
		while(own->p->X[0]>=start[0]+cell_size[0]*(ic[0]+1)){
			ic[0]++;
		}
		ic[0]+=ic_start[0];
		while(own->p->X[1]>=start[1]+cell_size[1]*(ic[1]+1)){
			ic[1]++;
		}
		ic[1]+=ic_start[1];
		tmp=own;
		if(own->next==NULL) {
			own=NULL;
		}else{
			own=own->next;
		}
		cells[local_index(ic)].insertParticle(tmp);
	}

	// If needed: output of created data
//	std::cout<<"Process "<<rank<<" received "<<num_part<<" particles\n";
//
//	if(num_part>0){
//		int n=0;
//		for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//			for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//				n+=cells[local_index(ic)].num_part;
//				std::cout<<"\nCell "<<local_index(ic)<<" of process "<<rank<<" contains "<<cells[local_index(ic)].num_part<<" particles\n";
//				for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
//					std::cout<<"Process "<<rank<<" received Particle at Position:\t"<<i->p->X[0]<<"/"<<i->p->X[1]<<"\n";
//				}
//			}
//		}
//	}
}

int SimProcess::index(int* p, int* n){
	int real_p[DIM];
	for(int d=0; d<DIM; d++){
		real_p[d]=p[d];
		if(real_p[d]<0){
			real_p[d] = n[d]+p[d];
		}else if(real_p[d]>=n[d]){
			real_p[d] = p[d]-n[d];
		}
	}

	#if 1==DIM
		return real_p[0];
	#endif

	#if 2==DIM
		return real_p[1]*n[0]+real_p[0];
	#endif

	#if 3==DIM
		return real_p[2]*n[0]*n[1]+real_p[1]*n[0]+real_p[0];
	#endif

	return 0;
}

int SimProcess::local_index(int* p){

	int loc_p[DIM];
	int loc_nc[DIM];
	for(int d=0; d<DIM; d++){
		loc_p[d]=p[d]-ic_start[d]+1;
		loc_nc[d]=local_nc[d]+2;
	}
//	std::cout<<"rank "<<rank<<" local_index returns "<<index(loc_p, loc_nc)<<"\n";
	return index(loc_p, loc_nc);
}

void SimProcess::status_output(std::string c){
	if(rank==0){
		std::cout<<c;
	}
}

void SimProcess::adding_to_pl(Cell* cells){
	int ic[DIM];
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			cells[local_index(ic)].adding_to_pl();
		}
	}
}

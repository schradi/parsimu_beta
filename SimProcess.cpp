/*
 * SImProcess.cpp
 *
 *  Created on: 21.04.2015
 *      Author: jonas
 */

#include "SimProcess.h"

void SimProcess::timeIntegration(real delta_t, real t_end, Cell* cells){
	// todo: Zeit stoppen, die zur Berechnung benötigt wird, zur Lastbalancierung interessant
	real t=0;
	max_V=cell_size[0]/delta_t;
	while (t<t_end){
//		compF();
//		compV();
//		compX();
//		moveParticles();
//		communicate();
		if(((int) (t/delta_t))%output_resolution==0){
			output(delta_t, cells);
		}
		t+=delta_t;
//		status_output("processing");
	}
}

void SimProcess::output(real delta_t, Cell* cells){
	MPI::COMM_WORLD.Barrier();
	int ic[DIM];
	if(rank==0) std::cout<<"\noutput:\n";
//	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//			for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
//				std::cout<<"Process "<<rank<<" received Particle at Position:\t"<<i->p->X[0]<<"/"<<i->p->X[1]<<"\n";
//			}
//		}
//	}

//	// Important: just particles stored in pl are considered

	MPI::Request request;
	MPI::Status status;

	int* send_pl_l = new int;
	*send_pl_l=get_num_p(ic_start, ic_stop, 0, 1, cells);
		/* direction of communication is not important,
		 * just needs to be different from other direction
		 */
	std::cout<<"Process "<<rank<<" contains "<<*send_pl_l<<" particles\n";

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
//					std::cout<<"Process "<<rank<<" received Particle at Position:\t"<<i->p->X[0]<<"/"<<i->p->X[1]<<"\n";
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

			std::cout<<"received length from p"<<cp<<"\t"<<recv_pl_l[cp]<<"\n";
			if(recv_pl_l[cp]!=0){
				std::cout<<"Particles found at process "<<cp<"\n";
				request=MPI::COMM_WORLD.Irecv(&recv_pl[pos], recv_pl_l[cp]*OUTP_SZE, MPI::DOUBLE, cp, global_np[0]*global_np[1]+cp);
				request.Wait(status);
			}
			std::cout<<"Received X of Particles:\n";
			for(int cc=1; cc<OUTP_SZE*recv_pl_l[cp]; cc+=OUTP_SZE){
				std::cout<<recv_pl[pos+cc]<<","<<recv_pl[pos+cc+1]<<"\n";
			}
			std::cout<<"\n";
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
		request=MPI::COMM_WORLD.Isend(send_pl_l, 1, MPI::INT, 0, 1);
		request.Wait(status);

		std::cout<<"length right now: "<<*send_pl_l<<"\n";

		// Sending Particles
		int ic[DIM];
		int pos=0;
		for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
			for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
				for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
//					std::cout<<"Process "<<rank<<" received Particle at Position:\t"<<i->p->X[0]<<"/"<<i->p->X[1]<<"\n";
					if (OUTP_SZE>0) send_pl[pos]=i->p->id;
					if (OUTP_SZE>1) send_pl[pos+1]=i->p->X[0];
					if (OUTP_SZE>2) send_pl[pos+2]=i->p->X[1];
					if (OUTP_SZE>3) send_pl[pos+3]=i->p->V[0];
					if (OUTP_SZE>4) send_pl[pos+4]=i->p->V[1];
					pos+=OUTP_SZE;

				}
			}
		}
		std::cout<<"P"<<rank<<" inserted till particle "<<pos/OUTP_SZE<<"\n";

//		code_pl(send_pl, ic_start, ic_stop, 0, 1, OUTP_SZE, cells);
		std::cout<<"Sended X[0] of Particles:\n";
		for(int cc=1; cc<*send_pl_l*OUTP_SZE; cc+=OUTP_SZE){
			std::cout<<send_pl[cc]<<"\n";
		}
		request=MPI::COMM_WORLD.Isend(send_pl, *send_pl_l*OUTP_SZE, MPI::DOUBLE, 0, global_np[0]*global_np[1]+rank);
		request.Wait(status);
	}
	delete send_pl_l;
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
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			c=&cells[local_index(ic)];
			for(ParticleList* pi=c->pl; pi->next!=NULL; pi=pi->next){
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
	int   c_id[DIM];
	real  X[DIM];
	bool  new_cell;
	ParticleList* pre;
	ParticleList* akt;
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			c=&cells[local_index(ic)];
			akt=c->pl;
			while(akt->next!=NULL){
				new_cell=false;
				for(int d=0; d<DIM; d++){
					new_ic[d]=ic[d];
					if(akt->p->X[d]<start[d]){
						new_ic[d]-=1;
						new_cell=true;
					}else if(akt->p->X[d]>=start[d]+cell_size[d]){
						new_ic[d]+=1;
						new_cell=true;
					}
				}
				if(new_cell){
					cells[local_index(new_ic)].insertParticle(akt->p);
					cells[local_index(ic)].deleteParticle(pre, akt);
				}
				pre=akt;
				akt=akt->next;
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

	if(rank==2){
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
	// For all surrounding Ghost-Cells: delete Pl
	int ic[DIM];
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			std::cout<<local_index(ic)<<"\t";
			cells[local_index(ic)].deletePl();
		}
		std::cout<<"\n";
	}
//	for(ic[1]=ic_start[1];ic[1]<=ic_stop[1]; ic[1]++){
//		ic[0]=ic_start[0]-1;
//		cells[local_index(ic)]->deletePl();
//		ic[0]=ic_stop[0]+1;
//		cells[local_index(ic)]->deletePl();
//	}


//	// For all cells: include adding in PL
//	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
//		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
//			cells[local_index(ic)]->adding_to_pl();
//		}
//	}

//	communicate(0, cells);
//	communicate(1, cells);
}

void SimProcess::communicate(int com_d, Cell* cells){
	int oth_d;
	if(com_d==0){
		oth_d=1;
	}else if(com_d==1){
		oth_d=0;
	}
	int ic[DIM];
//	std::vector<real> send_pl;
//	std::vector<real> recv_pl;
	int icr_start[DIM];
	int icr_stop[DIM];

	// Definition: 	right = ascending in com_d direction
	// 				left = descending in com_d direction

	/**
	 * communication to the right
	 */

	// Collect all pl within the cells of the right border and put them in an array of double
	icr_start[com_d]=ic_stop[com_d];
	icr_start[oth_d]=ic_start[oth_d]-1;
	icr_stop[com_d]=ic_stop[com_d]+1;
	icr_stop[oth_d]=ic_stop[oth_d]+1;

	int send_plr_l = get_num_p(icr_start, icr_stop, com_d, oth_d, cells);
	real send_plr[send_plr_l];

	int recv_pll_l;
	code_pl(send_plr, icr_start, icr_stop, com_d, oth_d, CD_P_SZE, cells);

	pl_l_send_recv(&send_plr_l, &recv_pll_l, com_d);

	real recv_pll[recv_pll_l];
	pl_send_recv(send_plr, send_plr_l, recv_pll, recv_pll_l, com_d);


	// decode array and sort into adding
	icr_start[com_d]=ic_start[com_d]-1;
	icr_start[oth_d]=ic_start[oth_d]-1;
	icr_stop[com_d]=ic_start[com_d];
	icr_stop[oth_d]=ic_stop[oth_d]+1;
	uncodePl(recv_pll, icr_start, icr_stop, com_d, oth_d, recv_pll_l, CD_P_SZE, cells);


	/**
	 * communication to the left
	 */

	// Collect all pl within the cells of the left border and put them in an array of double
	icr_start[com_d]=ic_start[com_d]-1;
	icr_start[oth_d]=ic_start[oth_d]-1;
	icr_stop[com_d]=ic_start[com_d];
	icr_stop[oth_d]=ic_stop[oth_d]+1;

	int send_pll_l = get_num_p(icr_start, icr_stop, com_d, oth_d, cells);
	real send_pll[send_pll_l];

	int recv_plr_l;
	code_pl(send_pll, icr_start, icr_stop, com_d, oth_d, CD_P_SZE, cells);

	pl_l_send_recv(&send_pll_l, &recv_plr_l, com_d);

	real recv_plr[recv_plr_l];
	pl_send_recv(send_pll, send_pll_l, recv_plr, recv_plr_l, com_d);


	// decode array and sort into adding
	icr_start[com_d]=ic_start[com_d]-1;
	icr_start[oth_d]=ic_start[oth_d]-1;
	icr_stop[com_d]=ic_start[com_d];
	icr_stop[oth_d]=ic_stop[oth_d]+1;
	uncodePl(recv_plr, icr_start, icr_stop, com_d, oth_d, recv_plr_l, CD_P_SZE, cells);

	// put all particles stored in adding in pl
	adding_to_pl(cells);
}

int SimProcess::get_num_p(int* icr_start, int* icr_stop, int com_d, int oth_d, Cell* cells){
	int ic[DIM];
	int num=0;
	if(rank==0) std::cout<<"--------------\n";
	for(ic[oth_d]=icr_start[oth_d]; ic[oth_d]<=icr_stop[oth_d]; ic[oth_d]++){
		for(ic[com_d]=icr_start[com_d]; ic[com_d]<=icr_stop[com_d]; ic[com_d]++){
			num+=cells[local_index(ic)].num_part;
		}
	}
	return num;
}

void SimProcess::pl_l_send_recv(int* send_pl_l, int* recv_pl_l, int com_d){
	MPI::Request request;
	MPI::Status status;

	// send length of the collected List
	request=MPI::COMM_WORLD.Isend(send_pl_l, 1, MPI::INT, neigh_upper[com_d], 1);

	// receive length of the received list
	request=MPI::COMM_WORLD.Irecv(recv_pl_l, 1, MPI::INT, neigh_lower[com_d], 1);
	request.Wait(status);
}

void SimProcess::pl_send_recv(real* send_pl, int send_pl_l, real* recv_pl, int recv_pl_l, int com_d){
	MPI::Request request;
	MPI::Status status;
	// send pl coded as array of double to the right
	request=MPI::COMM_WORLD.Isend(send_pl, send_pl_l, MPI::DOUBLE, neigh_upper[com_d], 2);

	// receive pl coded as array of double from the leftrsion from ‘int’ to ‘int*’ [-fpermissive]

	request=MPI::COMM_WORLD.Irecv(recv_pl, recv_pl_l, MPI::DOUBLE, neigh_lower[com_d], 2);

	// wait till communication is completed
	request.Wait(status);
}

void SimProcess::uncodePl(real* recv_pl, int* icr_start, int*icr_stop, int com_d, int oth_d, int length_recv, int size, Cell* cells){
	int pos=0;
	Particle* p;
	int ic[DIM];
	ic[oth_d]=icr_start[oth_d];
	ic[com_d]=icr_start[com_d];
	while(pos<length_recv){
		if(recv_pl[pos+2+com_d]<cells[local_index(ic)].start[com_d]+cell_size[com_d]){
			if(recv_pl[pos+3-com_d]<cells[local_index(ic)].start[oth_d]+cell_size[oth_d]){
				p=new Particle;
				cells[local_index(ic)].uncodePl(p, recv_pl, pos, size);
				pos+=size;
			}
		}else{
			if(ic[com_d]==ic_start[com_d]-1){
				ic[com_d]=ic_start[com_d];
			}else if(ic[com_d]==ic_start[com_d]){
				ic[com_d]=ic_start[com_d]-1;
				ic[oth_d]++;
			}
		}
	}
}

void SimProcess::code_pl(real* send_pl, int* icr_start, int* icr_stop, int com_d, int oth_d, int size, Cell* cells){
	int ic[DIM];
	int pos=0;
	for(ic[oth_d]=icr_start[oth_d]; ic[oth_d]<icr_stop[oth_d]; ic[oth_d]++){
		for(ic[com_d]=icr_start[com_d]; ic[com_d]<icr_stop[com_d]; ic[com_d]++){
			cells[local_index(ic)].code_pl(send_pl, pos, size);
			pos+=size;
		}
	}
}

void SimProcess::create_particles(ParticleList* new_pl, real* r_start, real* r_stop, real resolution){
	real pos[DIM];
	ParticleList* pl_iterate=new_pl;
	for(pos[1]=r_start[1]; pos[1]<r_stop[1]; pos[1]+=resolution){
		for(pos[0]=r_start[0]; pos[0]<r_stop[0]; pos[0]+=resolution){
			pl_iterate->p=new Particle;
			for(int d=0; d<DIM; d++){
				pl_iterate->p->X[d]=pos[d];
				pl_iterate->p->V[d]=0;
			}
			pl_iterate->p->V[1]=0.2;
			if(!(pos[1]>=r_stop[1]&&pos[0]>=r_stop[0])){
				pl_iterate->next=new ParticleList;
				pl_iterate=pl_iterate->next;
			}else{
				pl_iterate->next=NULL;
			}
		}
	}
}

void SimProcess::initData(Cell* cells){
	MPI::Request request;
	MPI::Status status;
	if(rank==0){
		// MASTER creates ParticleList containing all needed information
		real r_start[DIM]={0, 0};
		real r_stop[DIM]={10, global_size[1]};
		real resolution=1;
		ParticleList* to_insert=new ParticleList;
		create_particles(to_insert, r_start, r_stop, resolution);

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

		while(to_insert!=NULL){
			if(to_insert->p==NULL){
				to_insert=NULL;
			}else{
				proc[0]=0;proc[1]=0;
				while(to_insert->p->X[0]>=(proc[0]+1)*local_size[0]){
					proc[0]++;
				}
				while(to_insert->p->X[1]>=(proc[1]+1)*local_size[1]){
					proc[1]++;
				}
//				std::cout<<"Position:\t"<<proc[0]<<" "<<proc[1]<<"\n";
				proc_idx=index(proc, global_np);
//				std::cout<<"Processor:\t"<<proc_idx<<"\n";
//				std::cout<<"to_insert\t"<<to_insert<<"\n";
//				std::cout<<"proc_idx:\t"<<proc_idx<<"\n";


				// current particel is placed in the beginning of the list for each process
				tmp=sorted_Particles[proc_idx];
				sorted_Particles[proc_idx]=to_insert;
				to_insert=to_insert->next;
				sorted_Particles[proc_idx]->next=tmp;

				*pl_l[proc_idx]+=1;
			}
		}
//		for(int ip=0; ip<global_np[0]*global_np[1]; ip++){
//			std::cout<<"Process "<<ip<<" will contain "<<*pl_l[ip]<<" particles\n";
//		}

//		// Particles are send to the different Processes
//		// Sending number of Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			request=MPI::COMM_WORLD.Isend(pl_l[p], 1, MPI::INT, p, 1);
			request.Wait(status);
		}
//		// Sending Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			real cd_pl[*pl_l[p]*CD_P_SZE];
			ParticleList* tmp;
			tmp=sorted_Particles[p];
//			std::cout<<"pl_l["<<p<<"]\t"<<*pl_l[p]<<"\n";
			for(int i=0; i<*pl_l[p]*CD_P_SZE; i+=CD_P_SZE){
				cd_pl[i]=(real) tmp->p->id;
				if(CD_P_SZE>1) cd_pl[1+i]=tmp->p->X[0];
				if(CD_P_SZE>2) cd_pl[2+i]=tmp->p->X[1];
				if(CD_P_SZE>3) cd_pl[3+i]=tmp->p->V[0];
				if(CD_P_SZE>4) cd_pl[4+i]=tmp->p->V[1];
				if(CD_P_SZE>5) cd_pl[5+i]=tmp->p->m;
				tmp=tmp->next;
			}
			request=MPI::COMM_WORLD.Isend(cd_pl, *pl_l[p]*CD_P_SZE, MPI::DOUBLE, p, 2);
			request.Wait(status);
		}
		int ic[DIM];
		tmp=sorted_Particles[0];
//		std::cout<<"Anzahl pl in 0:\t"<<*pl_l[0]<<"\n";
		for(int pos=0; pos<*pl_l[0]; pos++){
			ic[0]=0;ic[1]=0;
			while(tmp->p->X[0]>=start[0]+cell_size[0]*ic[0]){
				ic[0]++;
			}
			ic[0]+=ic_start[0]-1;
			while(tmp->p->X[1]>=start[1]+cell_size[1]*ic[1]){
				ic[1]++;
			}
			ic[1]+=ic_start[1]-1;
			sorted_Particles[0]=sorted_Particles[0]->next;
			tmp->next=cells[local_index(ic)].pl;
			cells[local_index(ic)].pl=tmp;
			tmp=sorted_Particles[0];
			cells[local_index(ic)].num_part++;
			num_part++;
//			std::cout<<"particle sorted in cell ("<<ic[0]<<"/"<<ic[1]<<") = "<<local_index(ic)<<"\n";
		}
	}else{
		// waiting to receive particles
		int* new_recv;
		new_recv=new int;
		*new_recv=0;
		request=MPI::COMM_WORLD.Irecv(new_recv, 1, MPI::INT, 0, 1);
		request.Wait(status);
		num_part=*new_recv;
//		std::cout<<"Nr Particles in "<<rank<<":\t"<<num_part<<"\n";
		delete new_recv;

		// Receiving Particlelists
		real cod_pl[CD_P_SZE*num_part];
		request=MPI::COMM_WORLD.Irecv(cod_pl, CD_P_SZE*num_part, MPI::DOUBLE, 0, 2);
		request.Wait(status);
		std::cout<<"Sucessull request\n";
		int ic[DIM];
		Particle* p;
		for(int pos=0; pos<CD_P_SZE*num_part; pos+=CD_P_SZE){
			ic[0]=0; ic[1]=0;
			while(cod_pl[pos+1]>=start[0]+cell_size[0]*ic[0]){
				ic[0]++;
			}
			ic[0]+=ic_start[0]-1;
			while(cod_pl[pos+2]>=start[1]+cell_size[1]*ic[1]){
				ic[1]++;
			}
			ic[1]+=ic_start[1]-1;
			p=new Particle();
			if(CD_P_SZE>0) p->id=cod_pl[pos];
			if(CD_P_SZE>1) p->X[0]=cod_pl[pos+1];
			if(CD_P_SZE>2) p->X[1]=cod_pl[pos+2];
			if(CD_P_SZE>3) p->V[0]=cod_pl[pos+3];
			if(CD_P_SZE>4) p->V[1]=cod_pl[pos+4];
			if(CD_P_SZE>5) p->m=cod_pl[pos+5];
			cells[local_index(ic)].insertParticle(p);
//			std::cout<<"Received Particle at Position:\t"<<p->X[0]<<"/"<<p->X[1]<<"\n";
//
//			cells[local_index(ic)]
		}
//		std::cout<<"\n\n<><><><><><><><><><><>"<< cells[5].pl->p->X[0];
	}
	int ic[DIM];

//	std::cout<<"Process "<<rank<<" received "<<num_part<<" particles\n";
//	int n=0;
//	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//			n+=cells[local_index(ic)].num_part;
//			std::cout<<"cell"<<local_index(ic)<<" of process "<<rank<<" contains "<<cells[local_index(ic)].num_part<<"particles\n";
//
////			for(ParticleList* i=cells[local_index(ic)].pl; i!=NULL; i=i->next){
////				std::cout<<"Process "<<rank<<" received Particle at Position:\t"<<i->p->X[0]<<"/"<<i->p->X[1]<<"\n";
////			}
//		}
//	}
//	std::cout<<"Cell in "<<rank<<" received "<<n<<" particles\n";
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
	return index(loc_p, loc_nc);
}

void SimProcess::status_output(std::string c){
	if(rank==0){
		std::cout<<c;
	}
}

void SimProcess::adding_to_pl(Cell* cells){
	int ic[DIM];
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]; ic[0]++){
			cells[local_index(ic)].adding_to_pl();
		}
	}
}

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
	// Important: just particles stored in pl are considered
	int com_d=0;
	int oth_d=1;

	MPI::Request request;
	MPI::Status status;

	int* send_pl_l = new int;
	*send_pl_l=get_num_p(ic_start, ic_stop, com_d, oth_d, cells);

	if(rank==0){
		int recv_pl_l[global_np[0]*global_np[1]];
		int recv_pl_al=*send_pl_l;
		recv_pl_l[0]=*send_pl_l;
		int* new_recv;
		new_recv=new int;
		for (int cp=1; cp<global_np[0]*global_np[1];cp++){
			*new_recv=0;
			request=MPI::COMM_WORLD.Irecv(new_recv, 1, MPI::INT, cp, 1);
			request.Wait(status);
			recv_pl_l[cp]=*new_recv;
			recv_pl_al+=*new_recv;
		}

		real recv_pl[OUTP_SZE*recv_pl_al];
		code_pl(recv_pl, ic_start, ic_stop, com_d, oth_d, OUTP_SZE, cells);
		int pos=recv_pl_l[0];
		for (int cp=1; cp<global_np[0]*global_np[1];cp++){
			request=MPI::COMM_WORLD.Irecv(&recv_pl[pos], recv_pl_l[cp], MPI::DOUBLE, cp, 2);
			request.Wait(status);
			pos+=recv_pl_l[cp];
		}

		std::fstream file;
		int outp_nr=t/(delta_t*output_resolution);
		char cline_nr[32];
		sprintf(cline_nr, "%d", outp_nr);
		char outputfile_new[20]="data_neu/data.csv.";
		strcat(outputfile_new, cline_nr);
		file.open(outputfile_new, std::ios::out|std::ios::trunc);
		for(int i=0; i<recv_pl_al; i++){
			file<<recv_pl[0]<<","<<recv_pl[1]<<","<<recv_pl[2]<<std::endl;
		}
		file.close();

	}else{
		real send_pl[*send_pl_l];
		request=MPI::COMM_WORLD.Isend(send_pl_l, 1, MPI::INT, 0, 1);
		request.Wait(status);

		code_pl(send_pl, ic_start, ic_stop, com_d, oth_d, OUTP_SZE, cells);

		request=MPI::COMM_WORLD.Isend(send_pl, *send_pl_l, MPI::DOUBLE, 0, 2);
		request.Wait(status);
	}
}

void SimProcess::compF(Cell* cells){
	int ic[DIM];
	int nc[DIM];
	Cell* ci;
	Cell* cj;
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			// for each Cell
			ci=&cells[local_index(ic, global_nc)];
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
			ci=&cells[local_index(ic, global_nc)];

			for(ParticleList* pi=ci->pl; pi->next!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for (nc[1]=ic[1]-1;nc[1]<=ic[1]+1;nc[1]++){
					for (nc[0]=ic[0]-1;nc[0]<=ic[0]+1;nc[0]++){
						// for each neighboring cell
						cj=&cells[local_index(nc, global_nc)];
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
			c=&cells[local_index(ic, global_nc)];
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
			c=&cells[local_index(ic, global_nc)];
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
			c=&cells[local_index(ic, global_nc)];
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
					cells[local_index(new_ic, global_nc)].insertParticle(akt->p);
					cells[local_index(ic, global_nc)].deleteParticle(pre, akt);
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
	return local_index(id, local_nc);
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
		while(pos+cell_size[d]<local_size[d]){
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
		ic_stop[d]=ic_start[d]+local_nc[d];
	}

	// Nachbarn berechnen
	if(ip[0]-1>=0){
		neigh_lower[0]=ip[0]-1+global_np[0]*ip[1];
	}else{
		neigh_lower[0]=global_np[0]*(ip[1]+1);
	}
	if(ip[0]+1<global_np[0]){
		neigh_upper[0]=ip[0]+1+global_np[0]*ip[1];
	}else{
		neigh_upper[0]=global_np[0]*ip[1]+1;
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

//	status_output((char*)"ready to start");
}

void SimProcess::create_cells(Cell* cells){
	// Create Cells
	real p_cell_start[DIM];
	int ic[DIM];
	for (ic[1]=ic_start[1]-1; ic[1]<ic_stop[1]+1; ic[1]++){
		for (ic[0]=ic_start[0]-1; ic[0]<ic_stop[0]+1; ic[0]++){
			for (int d=0; d<DIM; d++){
				p_cell_start[d]=start[d]+ic[d]*r_cut;
			}
			cells[local_index(ic, global_nc)].set_params(p_cell_start, index(ic, global_nc), r_cut);
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
			std::cout<<local_index(ic, global_nc)<<"\t";
			cells[local_index(ic, global_nc)].deletePl();
		}
		std::cout<<"\n";
	}
//	for(ic[1]=ic_start[1];ic[1]<=ic_stop[1]; ic[1]++){
//		ic[0]=ic_start[0]-1;
//		cells[local_index(ic, global_nc)]->deletePl();
//		ic[0]=ic_stop[0]+1;
//		cells[local_index(ic, global_nc)]->deletePl();
//	}


//	// For all cells: include adding in PL
//	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
//		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
//			cells[local_index(ic, global_nc)]->adding_to_pl();
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
	for(ic[oth_d]=icr_start[oth_d]; ic[oth_d]<=icr_stop[oth_d]; ic[oth_d]++){
		for(ic[com_d]=icr_start[com_d]; ic[com_d]<=icr_stop[com_d]; ic[com_d]++){
			num+=cells[local_index(ic, global_nc)].num_part;
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
		if(recv_pl[pos+2+com_d]<cells[local_index(ic, global_nc)].start[com_d]+cell_size[com_d]){
			if(recv_pl[pos+3-com_d]<cells[local_index(ic, global_nc)].start[oth_d]+cell_size[oth_d]){
				p=new Particle;
				cells[local_index(ic, global_nc)].uncodePl(p, recv_pl, pos, size);
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
	for(ic[oth_d]=icr_start[oth_d]; ic[oth_d]<=icr_stop[oth_d]; ic[oth_d]++){
		for(ic[com_d]=icr_start[com_d]; ic[com_d]<=icr_stop[com_d]; ic[com_d]++){

			cells[local_index(ic, global_nc)].code_pl(send_pl, pos, size);
			pos+=size;
		}
	}
}

void SimProcess::initData(Cell* cells){
	// One Particle in each Cell-Center
	int ic[DIM];
	Particle* p;
	for (ic[1]=ic_start[1]; ic[1]<ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<ic_stop[0]; ic[0]++){
			p=cells[local_index(ic, global_nc)].pl->p;
			p=new Particle;
			p->m=1;
			p->id=index(ic, global_nc);
			for(int d=0; d<DIM; d++){
				p->V[d]=1;
				p->X[d]=start[d]+0.5*cell_size[d];
			}
		}
	}
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

int SimProcess::local_index(int* p, int* n){
	int loc_p[DIM];
	for(int d=0; d<DIM; d++){
		loc_p[d]=p[d]-ic_start[d]+1;
	}
	return index(loc_p, local_nc);
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
			cells[local_index(ic, global_nc)].adding_to_pl();
		}
	}
}

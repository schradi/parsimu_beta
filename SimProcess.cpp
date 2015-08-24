/*
 * SImProcess.cpp
 *
 *  Created on: 21.04.2015
 *      Author: jonas
 */

#include "SimProcess.h"

void SimProcess::timeIntegration(Cell* cells, real* p_map){
//	for(int p=0; p<np; p++)	{
//		if(rank==p){
//			std::cout<<"\nRank"<<rank<<"\n";
//			for(int d=0; d<DIM; d++){
//				std::cout<<"P"<<rank<<"cell_size["<<d<<"]="<< cell_size[d]<<"\n";
//				std::cout<<"P"<<rank<<"global_nc["<<d<<"]="<< global_nc[d]<<"\n";
//				std::cout<<"P"<<rank<<"global_size["<<d<<"]="<< global_size[d]<<"\n";
//				std::cout<<"P"<<rank<<"ic_start["<<d<<"]="<< ic_start[d]<<"\n";
//				std::cout<<"P"<<rank<<"ic_stop["<<d<<"]="<< ic_stop[d]<<"\n";
//				std::cout<<"P"<<rank<<"ip["<<d<<"]="<< ip[d]<<"\n";
//				std::cout<<"P"<<rank<<"local_nc["<<d<<"]="<< local_nc[d]<<"\n";
//				std::cout<<"P"<<rank<<"local_size["<<d<<"]="<< local_size[d]<<"\n";
//				std::cout<<"P"<<rank<<"neigh_lower["<<d<<"]="<< neigh_lower[d]<<"\n";
//				std::cout<<"P"<<rank<<"neigh_upper["<<d<<"]="<< neigh_upper[d]<<"\n";
//				std::cout<<"P"<<rank<<"global_np["<<d<<"]="<< global_np[d]<<"\n";
//				std::cout<<"P"<<rank<<"start["<<d<<"]="<< start[d]<<"\n";
//			}
//			std::cout<<"P"<<rank<<"lj_force_r_cut="<<lj_force_r_cut<<"\n";
//			std::cout<<"P"<<rank<<"delta_t="<< delta_t<<"\n";
//			std::cout<<"P"<<rank<<"t_end="<< t_end<<"\n";
//			std::cout<<"P"<<rank<<"global_num_part="<<global_num_part <<"\n";
//			std::cout<<"P"<<rank<<"num_part="<<num_part <<"\n";
//			std::cout<<"P"<<rank<<"rank="<<rank <<"\n";
//			std::cout<<"P"<<rank<<"r_cut="<<r_cut <<"\n";
//			std::cout<<"P"<<rank<<"max_V="<<max_V <<"\n";
//			std::cout<<"P"<<rank<<"sigma="<<sigma <<"\n";
//			std::cout<<"P"<<rank<<"epsilon="<< epsilon<<"\n";
//			std::cout<<"P"<<rank<<"Vvar="<<Vvar <<"\n";
//			std::cout<<"P"<<rank<<"output_resolution="<<output_resolution <<"\n";
//			std::cout<<"P"<<rank<<"t="<<t <<"\n";
//			std::cout<<"P"<<rank<<"np="<< np<<"\n";
//			std::cout<<"P"<<rank<<"output_folder="<<output_folder <<"\n";
//			std::cout<<"P"<<rank<<"log_time="<< log_time<<"\n";
//			std::cout<<"P"<<rank<<"log_energy="<< log_energy<<"\n";
//			std::cout<<"P"<<rank<<"log_positions="<< log_positions<<"\n";
//			std::cout<<"P"<<rank<<"log_velocity="<< log_velocity<<"\n";
//			std::cout<<"P"<<rank<<"log_id="<< log_id<<"\n";
//		}
//		MPI::COMM_WORLD.Barrier();
//	}
	t=delta_t;
	long int t_step_nr=1;
//	if(DOKU>=0) if(rank==0) std::cout<<"-----------------------\n- starting Simulation -\n-----------------------\n";
	compA(cells);
	while (t<t_end){
		clock_t t_start;
		if(rank==0) t_start=clock();
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
//		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compX\n";
		compX(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
//		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - moveParticles\n";
		moveParticles(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
//		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - communicate\n";
//		if(rank == 0) if((t_step_nr*100%(int)(t_end/delta_t))==0) std::cout<<"GHOST Part Vorher: "<<num_ghost_part<<"\n";
		communicate(cells);
//		if(rank == 0) if((t_step_nr*100%(int)(t_end/delta_t))==0) std::cout<<"GHOST Part Nachher: "<<num_ghost_part<<"\n";
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
//		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compA\n";
		compA(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
//		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compV\n";
		compV(cells);
		if(t_step_nr%output_resolution==0){
			if(DOKU>=2) std::cout<<"Pr "<<rank<<" - output\n";
			output(cells, (int) t_step_nr/output_resolution);
		}
		if(rank == 0) if((t_step_nr*100%(int)(t_end/delta_t))==0) std::cout<<"Process: "<<(int)((t/t_end)*100) + 1 <<"%\n";
		t+=delta_t;
		t_step_nr++;
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(rank==0 && log_time) timer_list->calc_avg_time("timeIntegration", t_start);
	}
}

void SimProcess::output(Cell* cells, int outp_nr){
	//compE(cells);
	if(log_positions){
		clock_t t_start;
		if(rank==0) t_start=clock();
		int ic[DIM];
//		for(ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//			for(ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//				for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
//					std::cout<<"P"<<rank<<" included Particle ("<<pi->p->X[0]<<"/"<<pi->p->X[1]<<") in Cell ("<<ic[0]<<"/"<<ic[1]<<")\n";
//				}
//			}
//		}
//		int c_tmp=0;
//		int ic[DIM];
//		for(ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//			for(ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//	//			num+=cells[local_index(ic)].num_part;
//				for(ParticleList* pi=cells[local_index(ic)].pl; pi pi=pi->next){
//					std::cout<<"P"<<rank<<" included Particle ("<<pi->p->X[0]<<"/"<<pi->p->X[1]<<") in Cell ("<<ic[0]<<"/"<<ic[1]<<")\n";
//				}
//			}
//		}

		// Important: just particles stored in pl are considered


//		if(DOKU>=3 )std::cout<<"Pr "<<rank<<" - output - begin\n";
		if(rank==0){

			// MASTER: Receives information from other processes and puts them in a file
			long int recv_pl_l[np]; // Number of particles of each process
			long int recv_pl_al=num_part;	// Number of particles at all

			// Receiving Number of Particles of each Process
			long int* new_recv;
			new_recv=new long int;
			recv_pl_l[0]=num_part;
			for (int cp=1; cp<MPI::COMM_WORLD.Get_size();cp++){
				*new_recv=0;
				MPI::COMM_WORLD.Recv(new_recv, 1, MPI::LONG, cp, 1);
				recv_pl_l[cp]=*new_recv;
				recv_pl_al+=*new_recv;
//				std::cout<<"received_length from "<<cp<<": "<<*new_recv<<"\n";
			}
			delete new_recv;
			// Receiving Particlelists
			real* recv_pl = new real[COM_SZE*recv_pl_al];

			code_range(recv_pl, ic_start, ic_stop, cells);
			long int pos=recv_pl_l[0]*COM_SZE;
			for (int cp=1; cp<np;cp++){
				MPI::COMM_WORLD.Recv(&recv_pl[pos], recv_pl_l[cp]*COM_SZE, MPI::DOUBLE, cp, 2);
				pos+=recv_pl_l[cp]*COM_SZE;
			}

			// Fileoutput
			std::fstream file;
			// Geting Name of the file
			char cline_nr[32];
			sprintf(cline_nr, "%d", outp_nr);
			char outputfile_new[20];
			strcpy(outputfile_new, output_folder);
			strcat(outputfile_new, "data.csv.");
			strcat(outputfile_new, cline_nr);
			file.open(outputfile_new, std::ios::out|std::ios::trunc);
		//		file<<"x coord, my coord, x velocity, y velocity\n";
			real corners[DIM];
			for(corners[1]=0; corners[1]<=global_size[1]; corners[1]+=global_size[1]){
				for(corners[0]=0; corners[0]<=global_size[0]; corners[0]+=global_size[0]){
					file<<corners[0]<<","<<corners[1];
					if(log_velocity) file<<",0,0";
					if(log_id) file<<",0";
					file<<"\n";
				}
			}
			for(pos=0; pos<recv_pl_al*COM_SZE; pos+=COM_SZE){
				/*
							 * 	code[0]=p->id;
								code[1]=p->m;
								code[2]=p->X[0];
								code[3]=p->V[0];
								code[4]=p->a_old[0];
								code[5]=p->X[1];
								code[6]=p->V[1];
								code[7]=p->a_old[1];
								code[8]=p->Epot;
								code[9]=p->Ekin;
							 */
				if(log_id) file<<recv_pl[pos]<<",";
				file<<recv_pl[pos+2]<<","<<recv_pl[pos+5];
				if(log_velocity) file<<","<<recv_pl[pos+3]<<","<<recv_pl[pos+6];
				file<<"\n";

			}
			file.close();
		}else{
			// Sending number of Particles
//			std::cout<<"Pr "<<rank<<" - sending nump: "<<num_part<<"\n";
			MPI::COMM_WORLD.Send(&num_part, 1, MPI::LONG, 0, 1);

			// Sending Particles
			real send_pl[num_part*COM_SZE];
			code_range(send_pl, ic_start, ic_stop, cells);
			MPI::COMM_WORLD.Send(send_pl, num_part*COM_SZE, MPI::DOUBLE, 0, 2);
		}
//		if(DOKU>=3) std::cout<<"Pr "<<rank<<" - output - end\n";
		if(rank==0 && log_time) timer_list->calc_avg_time("output", t_start);
	}
}

void SimProcess::compA(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	int ic[DIM];
	int nc[DIM];
	Cell* ci;
	Cell* cj;
	real space[DIM];
	Particle* p_own;
	Particle* p_oth;
	real F[DIM];
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			// for each Cell
			ci=&cells[local_index(ic)];
			for(ParticleList* pi=ci->pl; pi!=NULL; pi=pi->next){
				F[0]=0;
				F[1]=0;
				// for each Particle within this Cell
				for (nc[1]=ic[1]-1;nc[1]<=ic[1]+1;nc[1]++){
					for (nc[0]=ic[0]-1;nc[0]<=ic[0]+1;nc[0]++){
						// for each neighboring cell
						cj=&cells[local_index(nc)];
						for(ParticleList* pj=cj->pl; pj!=NULL; pj=pj->next ){
							if(pj!=pi){
								p_own=pi->p;
								p_oth=pj->p;
								space[0]=p_own->X[0]-p_oth->X[0];
								space[1]=p_own->X[1]-p_oth->X[1];
								force(space, F);
							}
						}
					}
				}
				for(int d=0; d<DIM; d++){
					pi->p->a_old[d]=pi->p->a[d];
					pi->p->a[d]=F[d]/pi->p->m;
				}
			}
		}
	}
	if(rank==0 && log_time) timer_list->calc_avg_time("compA", t_start);
}

void SimProcess::compV(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	Cell* c;
	int ic[DIM];
	real V_ges;
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			// for each Cell
			c=&cells[local_index(ic)];
			for(ParticleList* pi=c->pl; pi!=NULL; pi=pi->next){
				V_ges=0;
				// for each Particle within this Cell
				for(int d=0; d<DIM; d++){
					pi->p->V[d]+=0.5*(pi->p->a[d]+pi->p->a_old[d])*delta_t;
					V_ges=pow(pi->p->V[d],2);
				}
				pi->p->Ekin=0.5*pi->p->m*V_ges;
			}
		}
	}
	if(rank==0 && log_time) timer_list->calc_avg_time("compV", t_start);
}

void SimProcess::compX(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	Cell* c;
	int ic[DIM];
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			// for each Cell
			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
				// for each Particle within this Cell
				for(int d=0; d<DIM; d++){
					pi->p->X[d]+=pi->p->V[d]*delta_t+0.5*pi->p->a[d]*pow(delta_t,2);
				}
			}
		}
	}
	if(rank==0 && log_time) timer_list->calc_avg_time("compX", t_start);
}

void SimProcess::moveParticles(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	Cell* c;
	int   ic[DIM];
	int   new_ic[DIM];
	bool  new_cell;
	ParticleList* tmp;
	ParticleList* prev;
	ParticleList* akt;
	num_part=get_num_p(ic_start, ic_stop, cells);
//	if(num_part!=get_num_p(ic_start, ic_stop, cells)){
//		std::cout<<"P"<<rank<<"-FALSCH: movV\n";
//		while(t!=0){
//			t++;
//		}
//	}
	for (ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
		for (ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
			prev=NULL;
			// for each Cell
			c=&cells[local_index(ic)];
			akt=c->pl;
			while(akt!=NULL){
				new_cell=false;
				for(int d=0; d<DIM; d++){
					new_ic[d]=floor(akt->p->X[d]/cell_size[d]);
					if(ic[d]!=new_ic[d]){
						new_cell=true;
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
					if(new_ic[0]<ic_start[0] || new_ic[1]<ic_start[1] ||
							new_ic[0]>ic_stop[0] || new_ic[1]>ic_stop[1]){
						num_part--;
						num_ghost_part++;
					}
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
	num_part=get_num_p(ic_start, ic_stop, cells);
	if(rank==0 && log_time) timer_list->calc_avg_time("moveParticles", t_start);
}

int SimProcess::getLocalCellId(real* X){
	int id[DIM];
	id[0]=(int) ceil((X[0]/r_cut)-1);	// Abrunden
	id[1]=(int) ceil((X[1]/r_cut)-1);
	return local_index(id);
}

void SimProcess::devide_symetric(real* p_map){
	// muss davor bestimmt sein:
	// needs rank, np, global_size, cell_size, global_nc, global_np
	// sets: p_map

	if(rank==0){
//		std::cout<<"rank"<<rank<<"\n";
//		std::cout<<"np"<<np<<"\n";
//		std::cout<<"global_size[0]"<<global_size[0]<<"\n";
//		std::cout<<"global_size[1]"<<global_size[1]<<"\n";
//		std::cout<<"cell_size[0]"<<cell_size[0]<<"\n";
//		std::cout<<"cell_size[1]"<<cell_size[1]<<"\n";
//		std::cout<<"global_nc[0]"<<global_nc[0]<<"\n";
//		std::cout<<"global_nc[1]"<<global_nc[1]<<"\n";
//		std::cout<<"global_np[0]"<<global_np[0]<<"\n";
//		std::cout<<"global_np[1]"<<global_np[1]<<"\n";

		real np_size[DIM];
		real size[DIM];
		for(int d=0; d<DIM; d++){
			np_size[d]=global_size[d]/global_np[d];
		}
		int p_nr[DIM];
		for (p_nr[1]=0; p_nr[1]<global_np[1]; p_nr[1]++){
			for (p_nr[0]=0; p_nr[0]<global_np[0]; p_nr[0]++){
				int p=index(p_nr, global_np);
				for(int d=0; d<DIM; d++){
					p_map[p*DIM+d]=round(p_nr[d]*np_size[d]/cell_size[d])*cell_size[d];
//					std::cout<<"p_map["<<p*DIM+d<<"]="<<p_map[p*DIM+d]<<"\n";
//					size=p_map[p*DIM+d]+(np_size[d]/cell_size[d])*cell_size[d];
				}
			}
		}
//		for(int d=0; d<DIM; d++){
//			int i=0;
//			while(size[d]<global_size[d]){
//				p_map[p*DIM]
//				size+=cell_size[d];
//			}
//		}

//		for(int p=0; p<2*np; p+=2){
//			std::cout<<"map."<<p<<"=["<<p_map[p]<<","<<p_map[p+1]<<"]\n";
//		}
	}



}

void SimProcess::spread_global_info(){
	int nr_var=15;
	if(rank==0){
		// Send Variables
		// Wrap Variables
		real vars_send[nr_var];
		vars_send[0]=r_cut;
		vars_send[1]=global_size[0];
		vars_send[2]=global_size[1];
		vars_send[3]=t_end;
		vars_send[4]=delta_t;
		vars_send[5]=output_resolution;
		vars_send[6]=global_np[0];
		vars_send[7]=global_np[1];
		vars_send[8]=cell_size[0];
		vars_send[9]=cell_size[1];
		vars_send[10]=global_nc[0];
		vars_send[11]=global_nc[1];
		vars_send[12]=max_V;
		vars_send[13]=sigma;
		vars_send[14]=epsilon;
		for(int p=1; p<np; p++){
			MPI::COMM_WORLD.Send(vars_send, nr_var, MPI::DOUBLE, p, 1);
		}
	}else{
		// Receive the Variables
		real vars_recv[nr_var];
		MPI::COMM_WORLD.Recv(vars_recv, nr_var, MPI::DOUBLE, 0, 1);

		// Unwrap Variables
		r_cut			=vars_recv[0];
		global_size[0]	=vars_recv[1];
		global_size[1]	=vars_recv[2];
		t_end			=vars_recv[3];
		delta_t			=vars_recv[4];
		output_resolution=vars_recv[5];
		global_np[0]	=vars_recv[6];
		global_np[1]	=vars_recv[7];
		cell_size[0]	=vars_recv[8];
		cell_size[1]	=vars_recv[9];
		global_nc[0]	=vars_recv[10];
		global_nc[1]	=vars_recv[11];
		max_V			=vars_recv[12];
		sigma			=vars_recv[13];
		epsilon			=vars_recv[14];
	}
}

void SimProcess::spread_local_info(real* p_map){

	// needs: p_map, global_np, rank
	// sends start and local_size to each process
	// sets ic_start, ic_stop, local_nc,
	real recv[DIM*2];
	if(rank==0){
//		std::cout<<"\n";
		real var_send[np*DIM*2];
		int p_idx[DIM];
		int p=0;
		int pos=0;
		int n_way[DIM];
		n_way[0]=1;
		n_way[1]=global_np[0];
		for(p_idx[1]=0; p_idx[1]<global_np[1]; p_idx[1]++){
			for(p_idx[0]=0; p_idx[0]<global_np[0]; p_idx[0]++){
				// bestimmen von start und berechnen von local_size
				for(int d=0; d<DIM; d++){
					// start
					var_send[pos]=p_map[DIM*p+d];
					pos++;
					// local_size
					if(p_idx[d]<global_np[d]-1){
						var_send[pos]=p_map[DIM*(p+n_way[d])+d]-p_map[DIM*p+d];
					}else{
						var_send[pos]=global_size[d]-p_map[DIM*p+d];
					}
//					std::cout<<"P"<<p<<" - local_size["<<d<<"]="<<var_send[pos]<<"\n";
					pos++;
				}
				p++;
			}
		}
//		std::cout<<"var_send[3]="<<var_send[3]<<"\n";
		for(p=1; p<np; p++){
			//send var_send, length=4
//			std::cout<<"sending to "<<p<<": "<<var_send[p*4]<<", "<<": "<<var_send[p*4+1]<<", "<<": "<<var_send[p*4+2]<<", "<<": "<<var_send[p*4+3]<<"\n";
			MPI::COMM_WORLD.Send(&var_send[p*4], 4, MPI::DOUBLE, p, 1);
		}
		for(int d=0; d<DIM*2; d++){
			recv[d]=var_send[d];
		}
	}else{
		// receive recv, length=4
		MPI::COMM_WORLD.Recv(recv, 4, MPI::DOUBLE, 0, 1);
	}
	for(int d=0; d<DIM; d++){
		start[d]=recv[d*DIM];
		local_size[d]=recv[d*DIM+1];
		ic_start[d]=(int) start[d]/ (int) cell_size[d];
		local_nc[d]=(int) local_size[d]/(int) cell_size[d];
		ic_stop[d]=ic_start[d]+local_nc[d]-1;
	}

//	for(int p=0; p<np; p++){
//		MPI::COMM_WORLD.Barrier();
//		if(rank==p)std::cout<<" ";
//		MPI::COMM_WORLD.Barrier();
//		if(rank==p){
//			std::cout<<"\nProcess "<<rank<<": \n";
//			std::cout<<"ic_start: ";
//			for(int d=0; d<DIM; d++)std::cout<<ic_start[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"ic_stop: ";
//			for(int d=0; d<DIM; d++)std::cout<<ic_stop[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"ip: ";
//			for(int d=0; d<DIM; d++)std::cout<<ip[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"local_nc: ";
//			for(int d=0; d<DIM; d++)std::cout<<local_nc[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"local_size: ";
//			for(int d=0; d<DIM; d++)std::cout<<local_size[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"global_np: ";
//			for(int d=0; d<DIM; d++)std::cout<<global_np[d]<<",";
//			std::cout<<"\n";
//			std::cout<<"start - "<<local_index(ic_start)<<" : ";
//			for(int d=0; d<DIM; d++)std::cout<<start[d]<<",";
//			std::cout<<"\n";
//		}
//	}
}
void SimProcess::calculate_local_constants(){
	//	Local calculated Variables
	// needs; global_np, rank, cell_size, delta_t
	// sets: ip, neigh_lower, neigh_upper, max_V

	if(cell_size[0]<cell_size[1]) {
		max_V=cell_size[0]/delta_t;
	}else{
		max_V=cell_size[1]/delta_t;
	}


	ip[0]=0;
	ip[1]=0;
	while(ip[0]+global_np[0]*ip[1]!=rank){
		if(ip[0]+1<global_np[0]){
			ip[0]++;
		}else{
			ip[0]=0;
			ip[1]++;
		}
	}

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
//		std::cout<<"neigh_upper[1]="<<global_np[1]<<"\n";
	}else{
		neigh_upper[1]=ip[0];
	}
}

SimProcess::SimProcess(char* p_output_folder){
	timer_list=new TimerList();
	// Setting initial Variables
	num_part=0;
	num_ghost_part=0;
	np = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	log_time = false;
	log_energy = false;
	log_id=false;
	output_folder=p_output_folder;
	if(strcmp(output_folder, "")==0){
		folder_output=false;
	}
	t=0;

	int nr_var=19;
	if(rank==0){
		// Input from File
		char f_path[40];
		if(folder_output){
			strcpy(f_path, output_folder);
			strcat(f_path, "settings.dat");
		}else{
			strcpy(f_path, "settings.dat");
		}

		std::fstream file;
//		std::cout<<"open "<<f_path<<"\n";
		file.open(f_path, std::ios::in);
		char line[200];
		double testvar=0;
		int myint=0;
		r_cut=-1;
		global_size[0]=-1;
		t_end=-1;
		delta_t=-1;
		output_resolution=-1;
		sigma=-1;
		epsilon=-1;
		Vvar=-1;
		while(file >> line){
			if(!strcmp(line, "r_cut")){
				file >> line;
				r_cut = atof(line);
			}
			if(!strcmp(line, "global_size[0]")){
				file >> line;
				global_size[0] = atof(line);
			}
			if(!strcmp(line, "global_size[1]")){
				file >> line;
				global_size[1] = atof(line);
			}
			if(!strcmp(line, "t_end")){
				file >> line;
				t_end = atof(line);
			}
			if(!strcmp(line, "delta_t")){
				file >> line;
				delta_t = atof(line);
			}
			if(!strcmp(line, "output_resolution")){
				file >> line;
				output_resolution = atof(line);
			}
			if(!strcmp(line, "sigma")){
				file >> line;
				sigma = atof(line);
			}
			if(!strcmp(line, "epsilon")){
				file >> line;
				epsilon = atof(line);
			}
			if(!strcmp(line, "Vvar")){
				file >> line;
				Vvar = atof(line);
			}
		}
		if(r_cut==-1) std::cout<< "r_cut was not declared\n";
		if(global_size[0]==-1) std::cout<< "global_size was not declared\n";
		if(t_end==-1) std::cout<< "t_end was not declared\n";
		if(delta_t==-1) std::cout<< "delta_t was not declared\n";
		if(output_resolution==-1) std::cout<< "output_resolution was not declared\n";
		if(sigma==-1) std::cout<< "sigma was not declared\n";
		if(epsilon==-1) std::cout<< "epsilon was not declared\n";
		if(Vvar==-1) std::cout<< "Vvar was not declared\n";

		file.close();
		// Set important global variables
		// Set cell_size
		for(int d=0; d<DIM; d++){
			global_nc[d]=1;
			while((global_nc[d]+1)*r_cut<=global_size[d]){
				global_nc[d]++;
			}
		}

		int akt_split[DIM];
		real akt_psize[DIM];
		real akt_rating=0;
		real best_rating=0;
		for(akt_split[0] = 1; akt_split[0]<=np; akt_split[0]++){
			if(np%akt_split[0]==0){
				akt_psize[0]=global_size[0]/akt_split[0];
				akt_split[1]=np/akt_split[0];
				akt_psize[1]=global_size[1]/akt_split[1];
				akt_rating=akt_psize[1]/akt_psize[0];
				if(akt_rating>1) akt_rating=1/akt_rating;
				if(akt_rating>best_rating){
					best_rating=akt_rating;
					for(int d=0; d<DIM; d++){
						global_np[d]=akt_split[d];
					}
				}
			}
		}

		for(int d=0; d<DIM; d++){
			cell_size[d]=global_size[d]/global_nc[d];
			if(cell_size[d]>global_size[d]/global_np[d]){
				cell_size[d]=global_size[d]/global_np[d];
				r_cut=cell_size[d];
			}
		}
		if(rank==0) std::cout<<"spliting the field in "<<global_np[0]<<"x"<<global_np[1]<<" Processes\n";
	}
	spread_global_info();
	calculate_local_constants();
	lj_force_r_cut=lj_force(r_cut);

	if(rank==0  && folder_output){
		// Output in file
		char f_path[40];
		strcpy(f_path, output_folder);
		strcat(f_path, "settings.dat");
		std::fstream file;


		file.open(f_path, std::ios::out | std::ios::trunc);
		file<<"---Included settings---\n";
		file<<"r_cut\t"<<r_cut<<"\n";
		file<<"global_size[0]\t"<<global_size[0]<<"\n";
		file<<"global_size[1]\t"<<global_size[1]<<"\n";
		file<<"t_end\t"<<t_end<<"\n";
		file<<"delta_t\t"<<delta_t<<"\n";
		file<<"output_resolution\t"<<output_resolution<<"\n";
		file<<"sigma\t"<<sigma<<"\n";
		file<<"epsilon\t"<<epsilon<<"\n";
		file<<"Vvar\t"<<Vvar<<"\n";

		file<<"\n---Calculated Global Settings---\n";
		file<<"force(r_cut)\t"<<lj_force(r_cut)<<"\n";
		for(int d=0; d<DIM; d++){
			file<<"global_np["<<d<<"]\t"<<global_np[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"global_nc["<<d<<"]\t"<<global_nc[d]<<"\n";
		}
		file<<"max_V\t"<<max_V<<"\n";

//		file<<"\n---Calculated Local Settings of rank 0---\n";
//		for(int d=0; d<DIM; d++){
//			file<<"local_size["<<d<<"]\t"<<local_size[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"cell_size["<<d<<"]\t"<<cell_size[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"local_nc["<<d<<"]\t"<<local_nc[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"neigh_lower["<<d<<"]\t"<<neigh_lower[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"neigh_upper["<<d<<"]\t"<<neigh_upper[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"start["<<d<<"]\t"<<start[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"ic_start["<<d<<"]\t"<<ic_start[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"ic_stop["<<d<<"]\t"<<ic_stop[d]<<"\n";
//		}
//		for(int d=0; d<DIM; d++){
//			file<<"ip["<<d<<"]\t\t"<<ip[d]<<"\n";
//		}
		file.close();
	}
//	std::cout<<"Process "<<rank<<" - ready to start\n";
}

void SimProcess::create_cells(Cell* cells){
	real p_cell_start[DIM];
	int ic[DIM];
	for (ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for (ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
//			if(rank==0) std::cout<<"P"<<rank<<"-Cell ["<<ic[0]<<","<<ic[1]<<"], start=(";
			for (int d=0; d<DIM; d++){
				p_cell_start[d]=ic[d]*r_cut;
//				if(rank==0) std::cout<<p_cell_start[d]<<",";
			}
//			if(rank==0) std::cout<<")\n";
			cells[local_index(ic)].set_params(p_cell_start, index(ic, global_nc), cell_size);
		}
	}
}

void SimProcess::force(real* X, real* F){
	real r=0; /**< length*/
	for (int d=0;d<DIM;d++){
		r+=X[d]*X[d];
	}
	r=sqrt(r);

	//calculate force
	if(r<=r_cut){
		// Normalization of X
		for (int d=0; d<DIM; d++){
			X[d]=X[d]/r;
		}
		real f = lj_force(r)-lj_force_r_cut;
		for (int d=0; d<DIM; d++){
			X[d] = f*X[d];
		}
	}else{
		X[0]=0;
		X[1]=0;
	}
	F[0]+=X[0];
	F[1]+=X[1];
}

real SimProcess::lj_force(real r){
	real f1;
	real f2;
	f1=pow(sigma/r, 12);
	f2=pow(sigma/r, 6);
	return 24*epsilon/r*(2*f1-f2);
}

void SimProcess::communicate(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	int ic[DIM];
	// For all surrounding Ghost-Cells: delete Pl
	// For all cells: include adding in PL
	Cell* c;
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			c=&cells[local_index(ic)];
			if(ic[0]==ic_start[0]-1 || ic[0] == ic_stop[0]+1 ||
					ic[1]==ic_start[1]-1 || ic[1] == ic_stop[1]+1){
				c->deletePl();
			}
			c->adding_to_pl();
		}
	}
	num_ghost_part=0;
	for(int d=0; d<DIM; d++){
		if(global_np[d]==1){
			copy_border_cells(d, cells);
		}else{
			communicate(d, cells);
		}
	}
//	std::cout<<"P"<<rank<<" after communicate: "<<num_part<<"\n";
//	std::cout<<"P"<<rank<<" after communicate: Ghost"<<num_ghost_part<<"\n";
//	if(rank==0){
//		for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
//			for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
//				c=&cells[local_index(ic)];
//				std::cout<<"Cell ["<<ic[0]<<","<<ic[1]<<"] start=("<<c->start[0]<<","<<start[1]<<")\n";
//				for(ParticleList* pi=c->pl; pi!= NULL; pi=pi->next){
//					std::cout<<"Cell ["<<ic[0]<<","<<ic[1]<<"] has particle at ("<<pi->p->X[0]<<","<<pi->p->X[1]<<")\n";
//				}
//			}
//		}
//	}
//	if(rank==0 && log_time) timer_list->calc_avg_time("communicate", t_start);
}

void SimProcess::copy_border_cells(int com_d, Cell* cells){
//	std::cout<<"HHHHHHHHHHHHHH\nHHHHHHHHHHHHHH\nHHHHHHHHHHHHHH\nHHHHHHHHHHHHHH\n";
	int oth_d;
	if(com_d==0){
		oth_d=1;
	}else if(com_d==1){
		oth_d=0;
	}
	int pb_corr;
	int icr_lower_start[DIM], icr_upper_start[DIM];
	int icr_lower_stop[DIM], icr_upper_stop[DIM];

	icr_lower_start[com_d]=ic_start[com_d]-1;
	icr_lower_start[oth_d]=ic_start[oth_d]-1;
	icr_lower_stop[com_d]=ic_start[com_d];
	icr_lower_stop[oth_d]=ic_stop[oth_d]+1;

	icr_upper_start[com_d]=ic_stop[com_d];
	icr_upper_start[oth_d]=ic_start[oth_d]-1;
	icr_upper_stop[com_d]=ic_stop[com_d]+1;
	icr_upper_stop[oth_d]=ic_stop[oth_d]+1;

	long int pl_lower_length = get_num_p(icr_lower_start, icr_lower_stop, cells);
//	std::cout<<"copying pl_lower_length="<<pl_lower_length<<"\n";
	real pl_lower[pl_lower_length*COM_SZE];
	code_range(pl_lower, icr_lower_start, icr_lower_stop, cells);
	pb_corr=1;
	delete_pl(icr_lower_start, icr_lower_stop, cells);
	uncode_in_range(pl_lower, icr_upper_start, icr_upper_stop, pl_lower_length, cells, pb_corr, com_d);

	long int pl_upper_length = get_num_p(icr_upper_start, icr_upper_stop, cells);
//	std::cout<<"copying pl_upper_length="<<pl_upper_length<<"\n";
	real pl_upper[pl_upper_length*COM_SZE];
	code_range(pl_upper, icr_upper_start, icr_upper_stop, cells);
	pb_corr=-1;
	uncode_in_range(pl_upper, icr_lower_start, icr_lower_stop, pl_upper_length, cells, pb_corr, com_d);

}

void SimProcess::communicate(int com_d, Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	MPI::Request request;
	MPI::Status status;
	int oth_d;
	if(com_d==0){
		oth_d=1;
	}else if(com_d==1){
		oth_d=0;
	}
	int pb_corr;
	int icr_lower_start[DIM], icr_upper_start[DIM];
	int icr_lower_stop[DIM], icr_upper_stop[DIM];

	icr_lower_start[com_d]=ic_start[com_d]-1;
	icr_lower_start[oth_d]=ic_start[oth_d]-1;
	icr_lower_stop[com_d]=ic_start[com_d];
	icr_lower_stop[oth_d]=ic_stop[oth_d]+1;

	icr_upper_start[com_d]=ic_stop[com_d];
	icr_upper_start[oth_d]=ic_start[oth_d]-1;
	icr_upper_stop[com_d]=ic_stop[com_d]+1;
	icr_upper_stop[oth_d]=ic_stop[oth_d]+1;
//	for(int i=0; i<np; i++){
//		if(i==rank){
//			std::cout<<"Rank "<<rank<<" - Com_d: "<<com_d<<"\n";
//			std::cout<<"P"<<rank<<" icr_lower_start["<<com_d<<"]="<<icr_lower_start[com_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_lower_start["<<oth_d<<"]="<<icr_lower_start[oth_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_lower_stop["<<com_d<<"]="<<icr_lower_stop[com_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_lower_stop["<<oth_d<<"]="<<icr_lower_stop[oth_d]<<"\n";
//
//			std::cout<<"P"<<rank<<" icr_upper_start["<<com_d<<"]="<<icr_upper_start[com_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_upper_start["<<oth_d<<"]="<<icr_upper_start[oth_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_upper_stop["<<com_d<<"]="<<icr_upper_stop[com_d]<<"\n";
//			std::cout<<"P"<<rank<<" icr_upper_stop["<<oth_d<<"]="<<icr_upper_stop[oth_d]<<"\n";
//		}
//		MPI::COMM_WORLD.Barrier();
//	}

	int b=0;

	long int send_pl_lower_length = get_num_p(icr_lower_start, icr_lower_stop, cells);
//	if(rank==b) std::cout<<"send_pl_lower_length="<<send_pl_lower_length<<"\n";
	long int recv_pl_upper_length;
	real send_pl_lower[send_pl_lower_length*COM_SZE];

	code_range(send_pl_lower, icr_lower_start, icr_lower_stop, cells);
	MPI::COMM_WORLD.Isend(&send_pl_lower_length, 1, MPI::LONG, neigh_lower[com_d], 1);
	MPI::COMM_WORLD.Recv(&recv_pl_upper_length, 1, MPI::LONG, neigh_upper[com_d], 1);
	delete_pl(icr_lower_start, icr_lower_stop, cells);
//	if(rank==b) std::cout<<"recv_pl_upper_length="<<recv_pl_upper_length<<"\n";

	real recv_pl_upper[recv_pl_upper_length*COM_SZE];
	MPI::COMM_WORLD.Isend(&send_pl_lower, send_pl_lower_length*COM_SZE, MPI::DOUBLE, neigh_lower[com_d], 2);
	MPI::COMM_WORLD.Recv(&recv_pl_upper, recv_pl_upper_length*COM_SZE, MPI::DOUBLE, neigh_upper[com_d], 2);
	if(neigh_upper[com_d]<rank){
		pb_corr=1;
	}else{
		pb_corr=0;
	}
//	if(rank==b) std::cout<<"before uncode\n";
	uncode_in_range(recv_pl_upper, icr_upper_start, icr_upper_stop, recv_pl_upper_length, cells, pb_corr, com_d);
//	if(rank==b) std::cout<<"after uncode\n";


	// sending upper, receiving lower
	long int send_pl_upper_length = get_num_p(icr_upper_start, icr_upper_stop, cells);
//	if(rank==b) std::cout<<"send_pl_upper_length="<<send_pl_upper_length<<"\n";
	long int recv_pl_lower_length;
	real send_pl_upper[send_pl_upper_length*COM_SZE];
	code_range(send_pl_upper, icr_upper_start, icr_upper_stop, cells);
	MPI::COMM_WORLD.Isend(&send_pl_upper_length, 1, MPI::LONG, neigh_upper[com_d], 3);
	MPI::COMM_WORLD.Recv(&recv_pl_lower_length, 1, MPI::LONG, neigh_lower[com_d], 3);
//	if(rank==b) std::cout<<"recv_pl_lower_length="<<recv_pl_lower_length<<"\n";

	real recv_pl_lower[recv_pl_lower_length*COM_SZE];
	request = MPI::COMM_WORLD.Isend(&send_pl_upper, send_pl_upper_length*COM_SZE, MPI::DOUBLE, neigh_upper[com_d], 4);
	MPI::COMM_WORLD.Recv(&recv_pl_lower, recv_pl_lower_length*COM_SZE, MPI::DOUBLE, neigh_lower[com_d], 4);
	if(neigh_lower[com_d]>rank){
		pb_corr=-1;
	}else{
		pb_corr=0;
	}
	uncode_in_range(recv_pl_lower, icr_lower_start, icr_lower_stop, recv_pl_lower_length, cells, pb_corr, com_d);
	if(rank==0 && log_time) timer_list->calc_avg_time("communicate_1d", t_start);
}

void SimProcess::delete_pl(int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			if(ic[0]>=ic_start[0]&& ic[1]>=ic_start[1] &&
				ic[0]<=ic_stop[0] && ic[1]<=ic_stop[1]){
				num_part-=cells[local_index(ic)].num_part;
			}else{
				num_ghost_part-=cells[local_index(ic)].num_part;
			}
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
//			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
//				num++;
////				std::cout<<"P"<<rank<<" included Particle ("<<pi->p->X[0]<<"/"<<pi->p->X[1]<<") in Cell ("<<ic[0]<<"/"<<ic[1]<<")\n";
//			}
		}
	}
	return num;
}

void SimProcess::uncode_in_range(real* recv_pl, int* icr_start, int*icr_stop, long int length_recv, Cell* cells, int pb_corr, int com_d){
	// Eingabe: Codierte Partikel die von einem anderen Process geschickt wurden. Aufgrund der Periodischen Randbedingungen kann es deshalb sein, dass die Position der einzutragenden Position entspricht.
	// Ausgabe: Eingetragene Partikel in den Zellen mit richtiger Position
	// Partiklposition muss umgeschrieben werden wenn
	// bei einem Process:
	long int pos=0;
	Particle* p;
	int ic[DIM];
	int count=0;
//	if(rank==0) std::cout<<"P"<<rank<<" num_part="<<num_part<<"\n";
	while(pos<length_recv*COM_SZE){
		p=new Particle();
		uncode_p(&recv_pl[pos], p);
		// Periodic Boundaries
		if(pb_corr!=0 && com_d!=-1){
			p->X[com_d]+=pb_corr*global_size[com_d];
		}
		for(int d=0; d<DIM; d++){
			ic[d]=floor(p->X[d]/cell_size[d]);
//			std::cout<<"floor("<<p->X[d]<<"/"<<cell_size[d]<<")="<<ic[d]<<"\n";
		}
//		std::cout<<"P"<<rank<<" cell ["<<ic[0]<<","<<ic[1]<<"]="<<local_index(ic)<<"\n";
//		std::cout<<"size"<<cells[local_index(ic)].cell_size[1]<<"\n";
		cells[local_index(ic)].insertParticle(p);
//		if(rank==0) std::cout<<"p X("<<p->X[0]<<","<<p->X[1]<<") - insertet in Cell["<<ic[0]<<","<<ic[1]<<"]";
		if(p->X[0]>=start[0] && p->X[0]<start[0]+local_size[0] &&
			p->X[1]>=start[1] && p->X[1]<start[1]+local_size[1]){
			num_part++;
		}else{
//			if(rank==0)std::cout<<"-ghost";
			num_ghost_part++;
		}
//		if(rank==0)std::cout<<"\n";
		pos+=COM_SZE;
		count++;
	}
//	if(rank==0) std::cout<<"P"<<rank<<" AFTER num_part="<<num_part<<"\n";
//	std::cout<<"P"<<rank<<" counts "<<count<<" Particles\n";
//	for(ic[1]=ic_start[1]; ic[1]<=ic_stop[1]; ic[1]++){
//		for(ic[0]=ic_start[0]; ic[0]<=ic_stop[0]; ic[0]++){
//			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
//				std::cout<<"P"<<rank<<" included Particle ("<<pi->p->X[0]<<"/"<<pi->p->X[1]<<") in Cell ("<<ic[0]<<"/"<<ic[1]<<")\n";
//			}
//		}
//	}
}

void SimProcess::uncode_p(real* code, Particle* p){
	p->id=code[0];
	p->m=code[1];
	p->X[0]=code[2];
	p->V[0]=code[3];
	p->a_old[0]=code[4];
	p->X[1]=code[5];
	p->V[1]=code[6];
	p->a_old[1]=code[7];
	p->Epot=code[8];
	p->Ekin=code[9];
}

void SimProcess::code_range(real* send_pl, int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	long int pos=0;
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
				code_p(pi->p, &send_pl[pos]);
				pos+=COM_SZE;
			}
		}
	}
}

void SimProcess::code_p(Particle* p, real* code){
	code[0]=p->id;
	code[1]=p->m;
	code[2]=p->X[0];
	code[3]=p->V[0];
	code[4]=p->a_old[0];
	code[5]=p->X[1];
	code[6]=p->V[1];
	code[7]=p->a_old[1];
	code[8]=p->Epot;
	code[9]=p->Ekin;
}

void SimProcess::create_particle_cloud(ParticleList* new_pl, real* r_start, real* r_stop, real resolution, real* p_V){
	real pos[DIM];
	ParticleList* tmp;
	long int id_c=0;
	srand(time(NULL));
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
//				new_pl->p->V[d]=p_V[d]+(Vvar-(float) (rand()) / ((float) (RAND_MAX/(Vvar*2))));
				new_pl->p->V[d]=0.5;
			}
		}
	}
	std::cout<<"\nIncluded Particles\t"<<id_c<<"\n";
}

void SimProcess::create_particles(ParticleList* new_pl){
	if(num_part==0){
		real pos[DIM];
		ParticleList* tmp;
		long int id_c=0;
		char cfile[30];
		char in;


		std::cout<<"Would you like to use an existing file? - (y/n) - ";
		std::cin>>in;
		if(in=='y' || in=='Y'){
			std::cout<<"Would you like to use a Particle file stored in \"data/\"? - (y/n) - ";
			std::cin>>in;
			if(in=='y' || in=='Y'){
				char nr[4];
				std::cout<<"Which one would you like to use? - (number) - ";
				std::cin>>nr;
				strcpy(cfile, "data/data.csv.");
				strcat(cfile, nr);
			}else{
				std::cout<<"Insert path to particle file: ";
				std::cin>>cfile;
			}
			std::fstream file;
			file.open(cfile, std::ios::in);
			long int id;
			real m;
			real X[2]={0,0};
			real V[2]={0,0};
			real a[DIM];
			real a_old[DIM];
			std::string line;

			// Jump over Corners
			getline(file, line);
			getline(file, line);
			getline(file, line);
			getline(file, line);
			while(std::getline(file, line)){
				line+=',';
				size_t pos=0;
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					id=atoi(val.c_str());
				}
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					m=atof(val.c_str());
				}
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					X[0]=atof(val.c_str());
				}
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					X[1] = atof(val.c_str());
				}
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					V[0] = atof(val.c_str());
				}
				if ((pos = line.find(',')) != std::string::npos)
				{
					std::string val = line.substr(0, pos);
					line = line.substr(pos + 1);
					V[1] = atof(val.c_str());
				}

				id_c++;
				if(new_pl->p!=NULL){
					tmp=new_pl->next;
					new_pl->next=new ParticleList;
					new_pl->next->next=tmp;
					new_pl->next->p=new_pl->p;
					new_pl->p=NULL;
				}
				new_pl->p= new Particle;
				new_pl->p->id=id;
				new_pl->p->m=1;
				for(int d=0; d<DIM; d++){
					new_pl->p->X[d]=X[d];
					new_pl->p->V[d]=V[d];
					new_pl->p->a[d]=a[d];
					new_pl->p->a_old[d]=a_old[d];
				}
			}
			std::cout<<"Included Particles from File\t"<<id_c<<"\n";
		}else{
			in='y';
			real r_start[DIM];
			real r_stop[DIM];
			real res;
			real V[DIM];
			while(in=='y' || in=='Y'){
				std::cout<<"Would you like to create a new Particle Cloud? - ";
				std::cin>>in;
				if(in=='y' || in=='Y'){
					for(int d=0; d<DIM; d++){
						r_start[d]=-1;
						while(r_start[d]<0 || r_start[d]>global_size[d]){
							std::cout<<"Insert start["<<d<<"] - ";
							std::cin>>r_start[d];
						}
					}
					for(int d=0; d<DIM; d++){
						r_stop[d]=-1;
						while(r_stop[d]<=r_start[d] || r_stop[d]>global_size[d]){
							std::cout<<"Insert stop["<<d<<"] - ";
							std::cin>>r_stop[d];
						}
					}
					for(int d=0; d<DIM; d++){
						std::cout<<"Insert Velocity["<<d<<"] - ";
						std::cin>>V[d];
					}
					std::cout<<"Insert Resolution - ";
					std::cin>>res;
					create_particle_cloud(new_pl, r_start, r_stop, res, V);
				}
			}

		}
	}else{
		real r_start[DIM];
		real r_stop[DIM];
		real res;
		real V[DIM];
		r_stop[0]=global_size[0];
		r_stop[1]=global_size[1];
		V[0]=0;
		V[1]=0;
		res = global_size[0]/num_part;
		num_part=0;
		r_start[0]=0.5*res;
		r_start[1]=0.5*res;
		create_particle_cloud(new_pl, r_start, r_stop, res, V);
	}
}



void SimProcess::initData(Cell* cells, real* p_map){
	ParticleList* to_insert=new ParticleList;
	if(rank==0){
		create_particles(to_insert);
	}
	MPI::COMM_WORLD.Barrier();
	spread_particles(to_insert, cells, p_map);
}

void SimProcess::spread_particles(ParticleList* pl, Cell* cells, real* p_map){
	MPI::Request request;
	MPI::Status status;
	ParticleList* own=NULL;

	if(rank==0){
		// Sort Particles in different ParticleLists
		ParticleList* sorted_particles[np];
		long int *pl_l[np];
		for(int p=0; p<np; p++){
			pl_l[p]=new long int;
			*pl_l[p]=0;
		}
		sort_particles(pl, p_map, sorted_particles, pl_l);

		// Sending number of Particles
		for(int p=1; p<np; p++){
			MPI::COMM_WORLD.Send(pl_l[p], 1, MPI::LONG, p, 0);
		}

		// Sending Particles
		for(int p=1; p<np; p++){
			real cod_pl[*pl_l[p]*COM_SZE];
			ParticleList* tmp;
			tmp=sorted_particles[p];
			for(long int i=0; i<*pl_l[p]*COM_SZE; i+=COM_SZE){
				code_p(tmp->p, &cod_pl[i]);
				tmp=tmp->next;
			}
			MPI::COMM_WORLD.Send(cod_pl, *pl_l[p]*COM_SZE, MPI::DOUBLE, p, 1);
		}
		real cod_pl[*pl_l[0]*COM_SZE];
		long int i=0;
		ParticleList* pi=sorted_particles[0];
		while(i<*pl_l[0]){
			code_p(pi->p, &cod_pl[i*COM_SZE]);
			pi=pi->next;
			i++;
		}
		uncode_in_range(cod_pl, ic_start, ic_stop, *pl_l[0], cells);
	}else{
		// Receiving number of Particles
		long int recv_l;
		MPI::COMM_WORLD.Recv(&recv_l, 1, MPI::LONG, 0, 0);
		// Receiving Particles as Array of real
		real cod_pl[recv_l*COM_SZE];
		MPI::COMM_WORLD.Recv(cod_pl, COM_SZE*recv_l, MPI::DOUBLE, 0, 1);
		uncode_in_range(cod_pl, ic_start, ic_stop, recv_l, cells);
	}
}

void SimProcess::sort_particles(ParticleList* pl, real* p_map, ParticleList** sorted_particles, long int** pl_l){
	int t_ip[DIM];
	int p;
	ParticleList* tmp;
	while(pl!=NULL){
		t_ip[0]=global_np[0]-1;
		t_ip[1]=global_np[1]-1;
		real akt_border=p_map[t_ip[0]*2];
		while(pl->p->X[0]<p_map[t_ip[0]*2]){
			t_ip[0]--;
		}
		while(pl->p->X[1]<p_map[t_ip[1]*(global_np[0]*2)+1]){
			t_ip[1]--;
		}
		p=index(t_ip, global_np);
		*pl_l[p]+=1;
//		std::cout<<"Particle ("<<pl->p->X[0]<<"/"<<pl->p->X[1]<<") in Proc: "<<t_ip[0]<<"-"<<t_ip[1]<< "="<<p<< "\n";
		tmp=sorted_particles[p];
		sorted_particles[p]=pl;
		pl=pl->next;
		sorted_particles[p]->next=tmp;
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

int SimProcess::local_index(int* p){
	int loc_p[DIM];
	int loc_nc[DIM];
	for(int d=0; d<DIM; d++){
		loc_p[d]=p[d]-ic_start[d]+1;
		loc_nc[d]=local_nc[d]+2;
	}
	return index(loc_p, loc_nc);
}

void SimProcess::adding_to_pl(Cell* cells){
	int ic[DIM];
	for(ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for(ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			cells[local_index(ic)].adding_to_pl();
		}
	}
}

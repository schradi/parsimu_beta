/*
 * SImProcess.cpp
 *
 *  Created on: 21.04.2015
 *      Author: jonas
 */

#include "SimProcess.h"

void SimProcess::timeIntegration(Cell* cells){
	clock_t t_start;
	if(rank==0) t_start=clock();
	t=delta_t;
	long int t_step_nr=1;
	if(DOKU>=0) if(rank==0) std::cout<<"-----------------------\n- starting Simulation -\n-----------------------\n";
	compA(cells);
	while (t<t_end && !aborted){
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compX\n";
		compX(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - moveParticles\n";
		moveParticles(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - communicate\n";
		communicate(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compA\n";
		compA(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
		if(DOKU>=2) std::cout<<"Pr "<<rank<<" - compV\n";
		compV(cells);
		#if DEBUG
			MPI::COMM_WORLD.Barrier();
		#endif
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
	}
	if(rank==0) timerList->calc_avg_time("timeIntegration", t_start);
}

void SimProcess::output(Cell* cells, int outp_nr){
	clock_t t_start;
	if(rank==0) t_start=clock();
	int ic[DIM];

//	// Important: just particles stored in pl are considered
	num_part=get_num_p(ic_start, ic_stop, cells);


	if(DOKU>=3 )std::cout<<"Pr "<<rank<<" - output - begin\n";
	if(rank==0){
		clock_t t_debug;

		// MASTER: Receives information from other processes and puts them in a file
		int recv_pl_l[global_np[0]*global_np[1]]; // Number of particles of each process
		int recv_pl_al=num_part;	// Number of particles at all

		// Receiving Number of Particles of each Process
		int* new_recv;
		new_recv=new int;
		recv_pl_l[0]=num_part;
		for (int cp=1; cp<MPI::COMM_WORLD.Get_size();cp++){
			*new_recv=0;
			MPI::COMM_WORLD.Recv(new_recv, 1, MPI::INT, cp, 1);
			recv_pl_l[cp]=*new_recv;
			recv_pl_al+=*new_recv;
		}
		delete new_recv;
		// Receiving Particlelists
		real recv_pl[COM_SZE*recv_pl_al];

		code_range(recv_pl, ic_start, ic_stop, cells);
		int pos=recv_pl_l[0]*COM_SZE;
		for (int cp=1; cp<global_np[0]*global_np[1];cp++){
			MPI::COMM_WORLD.Recv(&recv_pl[pos], recv_pl_l[cp]*COM_SZE, MPI::DOUBLE, cp, 2);
			pos+=recv_pl_l[cp]*COM_SZE;
		}

		// Fileoutput
		std::fstream file;
		// Geting Name of the file
		char cline_nr[32];
		sprintf(cline_nr, "%d", outp_nr+1);
		char outputfile_new[20]="data/data.csv.";
		strcat(outputfile_new, cline_nr);
		file.open(outputfile_new, std::ios::out|std::ios::trunc);
//		file<<"x coord, my coord, x velocity, y velocity\n";
		file<<"0,0,0,0,0,0,0,0,0\n";
		file<<"0,0,"<<global_size[0]<<",0,0,0,0,0,0\n";
		file<<"0,0,0,0,0,"<<global_size[1]<<",0,0,0\n";
		file<<"0,0,"<<global_size[0]<<",0,0,"<<global_size[1]<<",0,0,0\n";
		real Ekin_ges=0;
		real Epot_ges=0;
		if(rank==0) t_debug=clock();
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
			file<<recv_pl[pos]<<","<<recv_pl[pos+2]<<","<<recv_pl[pos+3]<<","<<recv_pl[pos+5]<<","<<recv_pl[pos+6]<<"\n";
//			Epot_ges+=recv_pl[pos+8];
//			Ekin_ges+=recv_pl[pos+9];

		}
		file.close();
		if(rank==0) timerList->calc_avg_time("t_debug", t_debug);
//		file.open("data/energy.csv", std::ios::out|std::ios::app);
//		file<<t<<","<<Ekin_ges+Epot_ges<<","<<Epot_ges<<","<<Ekin_ges<<"\n";
//		file.close();

	}else{
		// Sending number of Particles
		real send_pl[num_part*COM_SZE];
		MPI::COMM_WORLD.Send(&num_part, 1, MPI::INT, 0, 1);

		// Sending Particles
		code_range(send_pl, ic_start, ic_stop, cells);
		MPI::COMM_WORLD.Send(send_pl, num_part*COM_SZE, MPI::DOUBLE, 0, 2);
	}
	if(DOKU>=3) std::cout<<"Pr "<<rank<<" - output - end\n";
	if(rank==0) timerList->calc_avg_time("output", t_start);
}

void SimProcess::compA(Cell* cells){
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

}

void SimProcess::compV(Cell* cells){
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

}

void SimProcess::compX(Cell* cells){
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

SimProcess::SimProcess(){
	timerList=new TimerList();
	// Setting initial Variables
	num_part=0;
	num_ghost_part=0;
	aborted=false;
	np = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	t=0;

	int nr_var=19;
	if(rank==0){
		// Input from File
		char f_path[40];
//		char cinput;
//		std::cout<<"Would you like to use settings stored in \"data/settings\"? (y/n) - ";
//		std::cin>>cinput;
//		if(cinput!='y' && cinput!='Y'){
//			std::cout<<"Insert path: ";
//			std::cin>>f_path;
//		}else{
//			strcpy(f_path, "data/settings");
//		}
		strcpy(f_path, "data/settings");

		std::fstream file;
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
			if(!strcmp(line, "global_size")){
				file >> line;
				global_size[0] = atof(line);
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
		for(int d=0; d<DIM; d++){
			global_np[d]= pow(np,(double) 1/DIM);
			local_size[d]=global_size[d]/global_np[d];
			local_nc[d]=1;
			while((local_nc[d]+1)*r_cut<=local_size[d]){
				local_nc[d]++;
			}
			cell_size[d]=local_size[d]/local_nc[d];
			global_nc[d]=local_nc[d]*global_np[d];
		}
		max_V=cell_size[0]/delta_t;
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
		vars_send[8]=local_size[0];
		vars_send[9]=local_size[1];
		vars_send[10]=cell_size[0];
		vars_send[11]=cell_size[1];
		vars_send[12]=local_nc[0];
		vars_send[13]=local_nc[1];
		vars_send[14]=global_nc[0];
		vars_send[15]=global_nc[1];
		vars_send[16]=max_V;
		vars_send[17]=sigma;
		vars_send[18]=epsilon;
		for(int p=1; p<np; p++){
			MPI::COMM_WORLD.Send(vars_send, nr_var, MPI::DOUBLE, p, 1);
		}
	}else{
		// Receive the Variables
		real vars_recv[nr_var];
		MPI::COMM_WORLD.Recv(vars_recv, nr_var, MPI::DOUBLE, 0, 1);
//		// Unwrap Variables
		r_cut			=vars_recv[0];
		global_size[0]	=vars_recv[1];
		global_size[1]	=vars_recv[2];
		t_end			=vars_recv[3];
		delta_t			=vars_recv[4];
		output_resolution=vars_recv[5];
		global_np[0]	=vars_recv[6];
		global_np[1]	=vars_recv[7];
		local_size[0]	=vars_recv[8];
		local_size[1]	=vars_recv[9];
		cell_size[0]	=vars_recv[10];
		cell_size[1]	=vars_recv[11];
		local_nc[0]		=vars_recv[12];
		local_nc[1]		=vars_recv[13];
		global_nc[0]	=vars_recv[14];
		global_nc[1]	=vars_recv[15];
		max_V			=vars_recv[16];
		sigma			=vars_recv[17];
		epsilon			=vars_recv[18];
	}
	//	Local calculated Variables
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
	for(int d=0; d<DIM; d++){
		start[d]=ip[d]*local_size[d];
		ic_start[d]=ip[d]*local_nc[d];
		ic_stop[d]=ic_start[d]+local_nc[d]-1;
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
	}else{
		neigh_upper[1]=ip[0];
	}
	lj_force_r_cut=lj_force(r_cut);

	if(rank==0){
		// Output in file
		std::fstream file;
		file.open("data/settings", std::ios::out | std::ios::trunc);
		file<<"---Included settings---\n";
		file<<"r_cut\t"<<r_cut<<"\n";
		file<<"global_size\t"<<global_size[0]<<"\n";
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

		file<<"\n---Calculated Local Settings of rank 0---\n";
		for(int d=0; d<DIM; d++){
			file<<"local_size["<<d<<"]\t"<<local_size[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"cell_size["<<d<<"]\t"<<cell_size[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"local_nc["<<d<<"]\t"<<local_nc[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"neigh_lower["<<d<<"]\t"<<neigh_lower[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"neigh_upper["<<d<<"]\t"<<neigh_upper[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"start["<<d<<"]\t"<<start[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"ic_start["<<d<<"]\t"<<ic_start[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"ic_stop["<<d<<"]\t"<<ic_stop[d]<<"\n";
		}
		for(int d=0; d<DIM; d++){
			file<<"ip["<<d<<"]\t\t"<<ip[d]<<"\n";
		}
		file.close();
	}
	std::cout<<"Process "<<rank<<" - ready to start\n";
}

void SimProcess::create_cells(Cell* cells){
	real p_cell_start[DIM];
	int ic[DIM];
	for (ic[1]=ic_start[1]-1; ic[1]<=ic_stop[1]+1; ic[1]++){
		for (ic[0]=ic_start[0]-1; ic[0]<=ic_stop[0]+1; ic[0]++){
			for (int d=0; d<DIM; d++){
				p_cell_start[d]=ic[d]*r_cut;
			}
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
	MPI::COMM_WORLD.Barrier();

	for(int d=0; d<DIM; d++){
		communicate(d, cells);
	}
	if(rank==0) timerList->calc_avg_time("communicate", t_start);
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

	int send_pl_lower_length = get_num_p(icr_lower_start, icr_lower_stop, cells);
	int recv_pl_upper_length;
	real send_pl_lower[send_pl_lower_length*COM_SZE];

	code_range(send_pl_lower, icr_lower_start, icr_lower_stop, cells);
	request = MPI::COMM_WORLD.Isend(&send_pl_lower_length, 1, MPI::INT, neigh_lower[com_d], 1);
	delete_pl(icr_lower_start, icr_lower_stop, cells);
	MPI::COMM_WORLD.Recv(&recv_pl_upper_length, 1, MPI::INT, neigh_upper[com_d], 1, status);
	real recv_pl_upper[recv_pl_upper_length*COM_SZE];
	request = MPI::COMM_WORLD.Isend(&send_pl_lower, send_pl_lower_length*COM_SZE, MPI::DOUBLE, neigh_lower[com_d], 2);
	MPI::COMM_WORLD.Recv(&recv_pl_upper, recv_pl_upper_length*COM_SZE, MPI::DOUBLE, neigh_upper[com_d], 2, status);
	if(neigh_upper[com_d]<rank){
		pb_corr=1;
	}else{
		pb_corr=0;
	}
	uncode_in_range(recv_pl_upper, icr_upper_start, icr_upper_stop, recv_pl_upper_length, cells, pb_corr, com_d);

	// sending upper, receiving lower
	int send_pl_upper_length = get_num_p(icr_upper_start, icr_upper_stop, cells);
	int recv_pl_lower_length;
	real send_pl_upper[send_pl_upper_length*COM_SZE];
	code_range(send_pl_upper, icr_upper_start, icr_upper_stop, cells);
	request = MPI::COMM_WORLD.Isend(&send_pl_upper_length, 1, MPI::INT, neigh_upper[com_d], 3);
	MPI::COMM_WORLD.Recv(&recv_pl_lower_length, 1, MPI::INT, neigh_lower[com_d], 3, status);
	request.Wait(status);

	real recv_pl_lower[recv_pl_lower_length*COM_SZE];
	request = MPI::COMM_WORLD.Isend(&send_pl_upper, send_pl_upper_length*COM_SZE, MPI::DOUBLE, neigh_upper[com_d], 4);
	MPI::COMM_WORLD.Recv(&recv_pl_lower, recv_pl_lower_length*COM_SZE, MPI::DOUBLE, neigh_lower[com_d], 4, status);
	request.Wait(status);

	if(neigh_lower[com_d]>rank){
		pb_corr=-1;
	}else{
		pb_corr=0;
	}
	uncode_in_range(recv_pl_lower, icr_lower_start, icr_lower_stop, recv_pl_lower_length, cells, pb_corr, com_d);
	if(rank==0) timerList->calc_avg_time("communicate_1d", t_start);
}

void SimProcess::delete_pl(int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			num_part-=cells[local_index(ic)].num_part;
			cells[local_index(ic)].deletePl();
		}
	}
}

int SimProcess::get_num_p(int* icr_start, int* icr_stop, Cell* cells){
	int ic[DIM];
	int num=0;
	for(ic[1]=icr_start[1]; ic[1]<=icr_stop[1]; ic[1]++){
		for(ic[0]=icr_start[0]; ic[0]<=icr_stop[0]; ic[0]++){
			for(ParticleList* pi=cells[local_index(ic)].pl; pi!=NULL; pi=pi->next){
				num++;
			}
		}
	}
	return num;
}

void SimProcess::uncode_in_range(real* recv_pl, int* icr_start, int*icr_stop, int length_recv, Cell* cells, int pb_corr, int com_d){
	// Eingabe: Codierte Partikel die von einem anderen Process geschickt wurden. Aufgrund der Periodischen Randbedingungen kann es deshalb sein, dass die Position der einzutragenden Position entspricht.
	// Ausgabe: Eingetragene Partikel in den Zellen mit richtiger Position
	// Partiklposition muss umgeschrieben werden wenn
	// bei einem Process:
	int pos=0;
	Particle* p;
	int ic[DIM];
	while(pos<length_recv*COM_SZE){
		p=new Particle();
		uncode_p(&recv_pl[pos], p);
		// Periodic Boundaries
		if(pb_corr!=0 && com_d!=-1){
			p->X[com_d]+=pb_corr*global_size[com_d];
		}
		ic[0]=icr_start[0];
		ic[1]=icr_start[1];
		for(int d=0; d<DIM; d++){
			while(p->X[d]>=cells[local_index(ic)].start[d]+cell_size[d]){
				ic[d]++;
			}
		}
//		std::cout<<"P"<<rank<<" cell ["<<ic[0]<<","<<ic[1]<<"]="<<local_index(ic)<<"\n";
//		std::cout<<"size"<<cells[local_index(ic)].cell_size[1]<<"\n";
		cells[local_index(ic)].insertParticle(p);
		if(p->X[0]>start[0] && p->X[0]<start[0]+local_size[0] &&
			p->X[1]>start[1] && p->X[1]<start[1]+local_size[1]){
			num_part++;
		}else{
			num_ghost_part++;
		}
		pos+=COM_SZE;
	}
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
	int pos=0;
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

void SimProcess::create_particles(ParticleList* new_pl, real* r_start, real* r_stop, real resolution, real* p_V){
	real pos[DIM];
	ParticleList* tmp;
	int id_c=0;
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
	std::cout<<"Included Particles\t"<<id_c<<"\n";


}

void SimProcess::insert_particles(ParticleList* new_pl){
	if(num_part==0){
		real pos[DIM];
		ParticleList* tmp;
		int id_c=0;
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
			int id;
			real m;
			real X[2]={0,0};
			real V[2]={0,0};
			real a[DIM];
			real a_old[DIM];
			std::string line;

			// Jump over first 5 Lines (Title + Corners)
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
					create_particles(new_pl, r_start, r_stop, res, V);
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
		create_particles(new_pl, r_start, r_stop, res, V);
	}
}

void SimProcess::initData(Cell* cells){
	MPI::Request request;
	MPI::Status status;
	ParticleList* own=NULL;
	int icr_start[DIM]={ic_start[0]-1, ic_start[1]-1};
	int icr_stop[DIM]={ic_stop[0]+1, ic_stop[1]+1};

	if(rank==0){
		ParticleList* to_insert=new ParticleList;
		insert_particles(to_insert);

		// Sort Particles in different ParticleLists
		ParticleList* tmp;
		ParticleList* sorted_Particles[global_np[0]*global_np[1]];
		int *pl_l[global_np[0]*global_np[1]];

		// Initialization
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

		// Sending number of Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			MPI::COMM_WORLD.Send(pl_l[p], 1, MPI::INT, p, p);
		}
		// Sending Particles
		for(int p=1; p<global_np[0]*global_np[1]; p++){
			real cod_pl[*pl_l[p]*COM_SZE];
			ParticleList* tmp;
			tmp=sorted_Particles[p];
			for(int i=0; i<*pl_l[p]*COM_SZE; i+=COM_SZE){
				code_p(tmp->p, &cod_pl[i]);
				tmp=tmp->next;
			}
			MPI::COMM_WORLD.Send(cod_pl, *pl_l[p]*COM_SZE, MPI::DOUBLE, p, global_nc[0]*global_nc[1]+p);
		}
		num_part=*pl_l[0];
		real cod_pl[num_part*COM_SZE];

		int i=0;

		for(ParticleList* pi=sorted_Particles[0]; pi!=NULL; pi=pi->next){
			code_p(pi->p, &cod_pl[i]);
			i+=COM_SZE;
		}
		uncode_in_range(cod_pl, icr_start, icr_stop, num_part, cells);
	}else{
		// Receiving number of Particles
		int recv_l;
		MPI::COMM_WORLD.Recv(&recv_l, 1, MPI::INT, 0, rank);
		num_part=recv_l;

		// Receiving Particles as Array of real
		real cod_pl[COM_SZE*num_part];
		MPI::COMM_WORLD.Recv(cod_pl, COM_SZE*num_part, MPI::DOUBLE, 0, global_nc[0]*global_nc[1]+rank);
		uncode_in_range(cod_pl, icr_start, icr_stop, num_part, cells);
	}
	// sort particles in Cells
//	ParticleList* tmp;
//	int ic[DIM];
//	while(own!=NULL){
//		ic[0]=0; ic[1]=0;
//		while(own->p->X[0]>=start[0]+cell_size[0]*(ic[0]+1)){
//			ic[0]++;
//		}
//		ic[0]+=ic_start[0];
//		while(own->p->X[1]>=start[1]+cell_size[1]*(ic[1]+1)){
//			ic[1]++;
//		}
//		ic[1]+=ic_start[1];
//		tmp=own;
//		if(own->next==NULL) {
//			own=NULL;
//		}else{
//			own=own->next;
//		}
//		cells[local_index(ic)].insertParticle(tmp);
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

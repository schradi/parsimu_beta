/*
 * timers.cpp
 *
 *  Created on: 08.07.2015
 *      Author: jonas
 */


#include "timers.h"

void Timer::calc_avg_time(clock_t start){
	clock_t curr=clock();
	timer=(timer*count+(curr-start))/(count+1);
	count++;
}

void Timer::calc_avg_time(clock_t start, clock_t curr){
	timer=(timer*count+(curr-start))/(count+1);
	count++;
}

void TimerList::calc_avg_time(std::string p_tag, clock_t start){
	bool found=false;
	clock_t curr=clock();
	TimerList* t=this;
	for(TimerList* ti=this; ti!=NULL; ti=ti->next){
		if(ti->t->tag==p_tag){
			found=true;
			ti->t->calc_avg_time(start, curr);
		}
		t=ti;
	}
	if(!found){
		t->next=new TimerList(p_tag);
	}
}

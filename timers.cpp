/*
 * timers.cpp
 *
 *  Created on: 08.07.2015
 *      Author: jonas
 */


#include "timers.h"

void Timer::calc_avg_time(time_t start){
	time_t curr;
	time(&curr);
	timer=(timer*count+(curr-start))/(count+1);
	count++;
}

void Timer::calc_avg_time(time_t start, time_t curr){
	timer=(timer*count+(curr-start))/(count+1);
	count++;
}

void TimerList::calc_avg_time(std::string p_tag, time_t start){
	bool found=false;
	time_t curr;
	time(&curr);
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

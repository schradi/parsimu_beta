/*
 * timers.h
 *
 *  Created on: 06.07.2015
 *      Author: jonas
 */

#ifndef TIMERS_H_
#define TIMERS_H_

#include <time.h>
#include <string>

class Timer{
public:
	clock_t timer;
	std::string tag;
	long int count;

	Timer (){
		timer=clock();
		count=0;
		tag="all";
	};
	Timer (std::string p_tag){
		timer=0;
		tag=p_tag;
		count=0;
	}
	void calc_avg_time(clock_t start);
	void calc_avg_time(clock_t start, clock_t curr);
};

class TimerList{
public:
	Timer* t;
	TimerList* next;

	TimerList(){
		t=new Timer();
		next=NULL;
	}
	TimerList(std::string p_tag){
		t=new Timer(p_tag);
		next=NULL;
	}
	void calc_avg_time(std::string p_tag, clock_t start);
};



#endif /* TIMERS_H_ */

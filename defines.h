/*
 * defines.h
 *
 *  Created on: 16.12.2014
 *      Author: jonas
 */

#ifndef DEFINES_H_
#define DEFINES_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

#define DOKU false
#define DET_DOKU false

#define DIM 2 /**< Dimension of the Simulation*/
#define CD_P_SZE 6
#define OUTP_SZE 3
#define NR_PROCS 4 /**< Number of Processes*/

//Used Datatype
typedef double real; /**< used datatype for all Values*/

#define sqr(aaa) ((aaa)*(aaa)) /**< definition of squaring*/

class Cell;
class Particle;
class Process;
class Simulation;

#endif /* DEFINES_H_ */
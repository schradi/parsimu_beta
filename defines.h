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

#define DET_DOKU false

#define DIM 2 /**< Dimension of the Simulation*/
#define COM_SZE 10
#define DOKU 1

#define DEBUG true

//Used Datatype
typedef double real; /**< used datatype for all Values*/

#define sqr(aaa) ((aaa)*(aaa)) /**< definition of squaring*/

class Cell;
class Particle;
class Process;
class Simulation;

#endif /* DEFINES_H_ */

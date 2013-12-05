#ifndef _ARGON_H_
#define _ARGON_H_
#include <math.h>

typedef struct
{
	double x, y, z;
} vect;

typedef struct
{
	vect
		*position,
		*momentum,
		*force;
	double
		potential,
		pressure,
		hamiltonian,
		temperature;
	
} state;

inline double len(vect position)
{
	return sqrt(position.x*position.x+position.y*position.y+position.z*position.z);
}

inline vect radius(vect first, vect second)
{
	first.x -= second.x;
	first.y -= second.y;
	first.z -= second.z;
	return first;
}

#endif

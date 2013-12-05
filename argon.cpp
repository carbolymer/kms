#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include "argon.hpp"
#include "ConfigFile.hpp"

using namespace std;

#define RANDOM_SIGN (drand48()>=0.5?1:-1)

/*
 * TODO: rozklady pedow przed symulacja i po symulacji - histogramy
 * TODO: trajektoria P-T dla kazdej z symulacji
 * TODO: wydrukowac sprawozdanie
 * TODO: dwa rozklady pedow: jednorodny i SB
 * */

int main(int argc, char *argv[])
{
	state systemState;
	double
		a,
		b0[3],
		b1[3],
		b2[3],
		T,
		m,
		f,
		L,
		tau,
		kb = 8.31e-3,
		R = 0.38e-9,
		epsilon = 1,
		initialEnergy,
		*pairPotentials,
		avTemperature = 0,
		avHamiltonian = 0,
		avPressure = 0,
		seconds = 0;		
	long int
		n,
		n3,
		n6,
		w = 0, i,j,k,nparticle,So,Sd,Sout,Sxyz;
	struct timeval timeBegin, timeStamp, newTimeStamp;

	string argument;

	bool isUniform = false;
		
	vect totalMomentum, *tMomentum;

	if(argc < 2)
	{
		cout << "No parameters file found!" << endl;
		return 1;
	}

	if(argc == 3)
	{
		argument = argv[2];
		if(argument.compare("uniform") == 0)
		{
			cout << "Using uniform momentum distribution." << endl;
			isUniform = true;
		}
	}

	ConfigFile parameters(argv[1]);

	n = parameters.getValueOfKey<int>("n");
	m = parameters.getValueOfKey<double>("m");
	epsilon = parameters.getValueOfKey<int>("epsilon");
	R = parameters.getValueOfKey<double>("R");
	f = parameters.getValueOfKey<double>("f");
	L = parameters.getValueOfKey<double>("L");
	a = parameters.getValueOfKey<double>("a");
	T = parameters.getValueOfKey<double>("T");
	tau = parameters.getValueOfKey<double>("tau");
	So = parameters.getValueOfKey<int>("So");
	Sd = parameters.getValueOfKey<int>("Sd");
	Sout = parameters.getValueOfKey<int>("Sout");
	Sxyz = parameters.getValueOfKey<int>("Sxyz");

	fstream crystalFile("avs.dat", fstream::out);
	fstream stateFile("state.dat", fstream::out);
	fstream initialFile("initial.dat", fstream::out);
	fstream finalFile("final.dat", fstream::out);

	cout << "Initial parameters:" << endl;
	cout << "n:\t\t" << n << endl; 
	cout << "m:\t\t" << m << endl; 
	cout << "epsilon:\t" << epsilon << endl; 
	cout << "R:\t\t" << R << endl; 
	cout << "f:\t\t" << f << endl; 
	cout << "L:\t\t" << L << endl; 
	cout << "a:\t\t" << a << endl; 
	cout << "T_0:\t\t" << T << endl; 
	cout << "tau:\t\t" << tau << endl;
	cout << "S_o:\t\t" << So << endl; 
	cout << "S_d:\t\t" << Sd << endl;
	cout << "S_out:\t\t" << Sout << endl; 
	cout << "S_xyz:\t\t" << Sxyz << endl << endl;

 	srand48(time(NULL));

	n3 = n*n*n;
	n6 = n3*n3;
	
	initialEnergy = -kb*T/2.;

	b0[0]=a;
	b0[1]=0.0;
	b0[2]=0.0;
	
	b1[0]=a*0.5;
	b1[1]=a*sqrt(3.)/2.;
	b1[2]=0.0;
	
	b2[0]=a*0.5;
	b2[1]=a*sqrt(3.)/6.;
	b2[2]=a*sqrt(2./3.);

	tMomentum = new vect[n3];
	
	systemState.position = new vect[n3];
	systemState.momentum = new vect[n3];
	systemState.force = new vect[n3];
	
	totalMomentum.x = 0;
	totalMomentum.y = 0;
	totalMomentum.z = 0;


	// initial conditions
	for(i=0; i<n; ++i)
		for(j=0; j<n; ++j)
			for(k=0; k<n; ++k)
			{
				systemState.position[w].x = (i-(n-1)*0.5)*b0[0]+(j-(n-1)*0.5)*b1[0]+(k-(n-1)*0.5)*b2[0];
				systemState.position[w].y = (i-(n-1)*0.5)*b0[1]+(j-(n-1)*0.5)*b1[1]+(k-(n-1)*0.5)*b2[1];
				systemState.position[w].z = (i-(n-1)*0.5)*b0[2]+(j-(n-1)*0.5)*b1[2]+(k-(n-1)*0.5)*b2[2];

				if(isUniform)
				{
					totalMomentum.x += systemState.momentum[w].x = sqrt(-2.*m*initialEnergy)*RANDOM_SIGN*drand48()*2.;
					totalMomentum.y += systemState.momentum[w].y = sqrt(-2.*m*initialEnergy)*RANDOM_SIGN*drand48()*2.;
					totalMomentum.z += systemState.momentum[w].z = sqrt(-2.*m*initialEnergy)*RANDOM_SIGN*drand48()*2.;
				}
				else
				{
					totalMomentum.x += systemState.momentum[w].x = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
					totalMomentum.y += systemState.momentum[w].y = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
					totalMomentum.z += systemState.momentum[w].z = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
				}

				++w;
			}
			
	totalMomentum.x /= n3;
	totalMomentum.y /= n3;
	totalMomentum.z /= n3;

	w = 0;
	for(i=0; i<n; ++i)
		for(j=0; j<n; ++j)
			for(k=0; k<n; ++k)	
			{
				systemState.momentum[w].x -= totalMomentum.x;
				systemState.momentum[w].y -= totalMomentum.y;
				systemState.momentum[w].z -= totalMomentum.z;
			
				++w;
			}

	//
	// **
	//		

	gettimeofday(&timeBegin,NULL);
	for(int step = 0; step < So + Sd ; ++step)
	{
		systemState.pressure = 0;
		systemState.potential = 0;
		systemState.temperature = 0;
		systemState.hamiltonian = 0;
					
		// calculations for pairs
		// i=j
		for(i = 0; i < n3; ++i)
		{
			systemState.force[i].x = 0;
			systemState.force[i].y = 0;
			systemState.force[i].z = 0;
		}
		// i > j
		for(i = 0; i < n3; ++i)
			for(j = 0; j < i; ++j)
			{
				vect vRadius = radius(systemState.position[i], systemState.position[j]);
				double distance = len(vRadius);
				double coefficient = 12.*epsilon*(pow(R/distance,12)-pow(R/distance,6))/distance/distance;
				
				systemState.potential += epsilon*(pow(R/distance,12)-2.*pow(R/distance,6));
				
				systemState.force[i].x += coefficient*vRadius.x;
				systemState.force[i].y += coefficient*vRadius.y;
				systemState.force[i].z += coefficient*vRadius.z;

				// j > i
				systemState.force[j].x += -coefficient*vRadius.x;
				systemState.force[j].y += -coefficient*vRadius.y;
				systemState.force[j].z += -coefficient*vRadius.z;
			}
			
		// calculation for particles
		for(i = 0; i < n3; ++i)
		{
			double distance = len(systemState.position[i]);
			double coefficient;
			if(distance >= L)
			{
				coefficient = f*(L-distance)/distance;

				systemState.force[i].x += coefficient*systemState.position[i].x;
				systemState.force[i].y += coefficient*systemState.position[i].y;
				systemState.force[i].z += coefficient*systemState.position[i].z;
				
				systemState.potential += f*pow(distance-L,2)/2.;
				systemState.pressure += len(systemState.force[i]);
			}		
		}
		systemState.pressure /= 4.*M_PI*L*L;
		
		// motion equations
		for(i = 0; i < n3; ++i)
		{
			double energy = 0;
			if(step > 0)
			{
				// 18c
				systemState.momentum[i].x = tMomentum[i].x + systemState.force[i].x*tau/2.;
				systemState.momentum[i].y = tMomentum[i].y + systemState.force[i].y*tau/2.;
				systemState.momentum[i].z = tMomentum[i].z + systemState.force[i].z*tau/2.;
			}

			energy = pow(len(systemState.momentum[i]),2)/2./m;
			systemState.temperature += 2./3./n3/kb * energy;
			systemState.hamiltonian += energy;

			// 18a
			tMomentum[i].x = systemState.momentum[i].x + systemState.force[i].x*tau/2.;
			tMomentum[i].y = systemState.momentum[i].y + systemState.force[i].y*tau/2.;
			tMomentum[i].z = systemState.momentum[i].z + systemState.force[i].z*tau/2.;

			// 18b
			systemState.position[i].x = systemState.position[i].x + systemState.momentum[i].x*tau/m;
			systemState.position[i].y = systemState.position[i].y + systemState.momentum[i].y*tau/m;
			systemState.position[i].z = systemState.position[i].z + systemState.momentum[i].z*tau/m;
		
			if(step % Sxyz == 0 && step > So)
			{
				crystalFile << systemState.position[i].x <<" "<< systemState.position[i].y <<" "<< systemState.position[i].z << " ";
				crystalFile << systemState.momentum[i].x <<" "<< systemState.momentum[i].y <<" "<< systemState.momentum[i].z << endl;
			}
			if(step == 0)
			{
				initialFile << systemState.position[i].x <<" "<< systemState.position[i].y <<" "<< systemState.position[i].z << " ";
				initialFile << systemState.momentum[i].x <<" "<< systemState.momentum[i].y <<" "<< systemState.momentum[i].z << endl;
			}
			else if(step == (So+Sd-1))
			{
				finalFile << systemState.position[i].x <<" "<< systemState.position[i].y <<" "<< systemState.position[i].z << " ";
				finalFile << systemState.momentum[i].x <<" "<< systemState.momentum[i].y <<" "<< systemState.momentum[i].z << endl;
			}
		}
		systemState.hamiltonian += systemState.potential;

		if(step > So)
		{
			avTemperature += systemState.temperature;
			avHamiltonian += systemState.hamiltonian;
			avPressure += systemState.pressure;
		}

		if(step % Sout == 0 && step > So)
		{
			stateFile << tau*step << " " << systemState.hamiltonian << " " << systemState.potential << " " << systemState.pressure << " " << systemState.temperature << endl;
			crystalFile << endl << endl;
		}

		gettimeofday(&timeStamp,NULL);
		if(step % 100 == 0 && step != 0)
		{
			seconds = (timeStamp.tv_sec - timeBegin.tv_sec + (timeStamp.tv_usec - timeBegin.tv_usec)/1e6f)/step*(So+Sd-step);
			cout << "Progress: \t" << step << "\t/ " << So + Sd << "\t\tElapsed time:\t" << (timeStamp.tv_sec - timeBegin.tv_sec) << " s";
			if(seconds > 0)
			{
				cout  << "\tRemaining time: " << seconds  << " s\t" << "= " << seconds/3600. << " h";
				if(seconds > 86400)
					cout << "\t= " << seconds/86400. << " d";
			}
		   	cout << endl;
		}
		if(step == 0)
			initialFile.close();
		if(step == (So+Sd-1))
			finalFile.close();
	}

	avTemperature /= Sd;
	avHamiltonian /= Sd;
	avPressure /= Sd;
	
	cout << "Simulation complete. " << endl;
	cout << "average H:\t" << avHamiltonian << " kJ/mol" << endl;
	cout << "average T:\t" << avTemperature << " K" << endl;
	cout << "average P:\t" << avPressure*16.6 << " atm" << endl;

	crystalFile.close();

	return 0;
}

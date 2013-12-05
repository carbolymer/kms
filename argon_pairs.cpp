#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "argon.hpp"
#include "ConfigFile.hpp"

using namespace std;

#define RANDOM_SIGN (drand48()>=0.5?1:-1)

/*
	polozenia i pedy do jednego pliku:
		x	y	z	px	py	pz
	kazdy krok czasowy oddzielic \n\n
	w drugim pliku zapisywac czas, potencjal i energie kinetyczna
	
	http://omega.if.pw.edu.pl/~mars/labKMS/
	
	dorobic parser plikow wejsciowych, sfera w wizualizacji

*/
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
		avPressure = 0;
	long int
		n,
		n3,
		n6,
		w = 0, i,j,k,nparticle,So,Sd,Sout,Sxyz;
		
	vect *pairDistances, *pairForces, totalMomentum, *tMomentum;

	if(argc < 1)
	{
		cout << "No parameters file found!" << endl;
		return 1;
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

	// fstream parametersFile;
	// if(argc > 1)
	// 	parametersFile.open(argv[1], fstream::in);
	// else
	// 	parametersFile.open("argon.in", fstream::in);

	// parametersFile >> n;		// 1
	// parametersFile >> m;		// 2
	// parametersFile >> epsilon;	// 3
	// parametersFile >> R;		// 4
	// parametersFile >> f;		// 5
	// parametersFile >> L;		// 6
	// parametersFile >> a;		// 7
	// parametersFile >> T;		// 8
	// parametersFile >> tau;		// 9
	// parametersFile >> So;		// 10
	// parametersFile >> Sd;		// 11
	// parametersFile >> Sout;		// 12
	// parametersFile >> Sxyz;		// 13

	// parametersFile.close();



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

	// vect *dupa = new vect[n6/n];
	// cout << "moge zaalokowac :" << sizeof(vect)*n6/1024./1024./1024.*4. << "GB" << endl;
	// return 1;
	pairDistances = new vect[n6]; // i+j*n3
	pairForces = new vect[n6];
	tMomentum = new vect[n3];
	pairPotentials = new double[n6];
	
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

				totalMomentum.x += systemState.momentum[w].x = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
				totalMomentum.y += systemState.momentum[w].y = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
				totalMomentum.z += systemState.momentum[w].z = sqrt(2.*m*initialEnergy*log(drand48()))*RANDOM_SIGN;
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
				// at the beginning system is not moving
				systemState.momentum[w].x -= totalMomentum.x;
				systemState.momentum[w].y -= totalMomentum.y;
				systemState.momentum[w].z -= totalMomentum.z;
			
				// crystalFile << systemState.position[w].x <<" "<< systemState.position[w].y <<" "<< systemState.position[w].z << " ";
				// crystalFile << systemState.momentum[w].x <<" "<< systemState.momentum[w].y <<" "<< systemState.momentum[w].z << endl;
				++w;
			}
	
	// crystalFile << endl;

	//
	// **
	//		

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
			nparticle = i+i*n3;
			pairDistances[nparticle].x = 0;
			pairDistances[nparticle].y = 0;
			pairDistances[nparticle].z = 0;
			
			pairForces[nparticle].x = 0;
			pairForces[nparticle].y = 0;
			pairForces[nparticle].z = 0;
			
			systemState.potential += pairPotentials[nparticle] = 0;

			systemState.force[i].x = 0;
			systemState.force[i].y = 0;
			systemState.force[i].z = 0;
		}
		// i > j
		for(i = 0; i < n3; ++i)
			for(j = 0; j < i; ++j)
			{
				int npair = i + j*n3;
				pairDistances[npair].x = systemState.position[i].x - systemState.position[j].x;
				pairDistances[npair].y = systemState.position[i].y - systemState.position[j].y;			
				pairDistances[npair].z = systemState.position[i].z - systemState.position[j].z;

				double distance = len(pairDistances[npair]);
				double coefficient = 12.*epsilon*(pow(R/distance,12)-pow(R/distance,6))/distance/distance;
		
				pairForces[npair].x = coefficient*pairDistances[npair].x;
				pairForces[npair].y = coefficient*pairDistances[npair].y;
				pairForces[npair].z = coefficient*pairDistances[npair].z;
				
				pairPotentials[npair] = epsilon*(pow(R/distance,12)-2.*pow(R/distance,6));
				
				systemState.force[i].x += pairForces[npair].x;
				systemState.force[i].y += pairForces[npair].y;
				systemState.force[i].z += pairForces[npair].z;
			}
		
		// i < j
		for(i = 0; i < n3; ++i)
			for(j = i+1; j < n3; ++j)
			{
				int npair = i + j*n3;
				int npair2 = j+i*n3;

				pairDistances[npair].x = -pairDistances[npair2].x;
				pairDistances[npair].y = -pairDistances[npair2].y;
				pairDistances[npair].z = -pairDistances[npair2].z;

				pairForces[npair].x = -pairForces[npair2].x;
				pairForces[npair].y = -pairForces[npair2].y;
				pairForces[npair].z = -pairForces[npair2].z;
				
				systemState.potential += pairPotentials[npair] = pairPotentials[npair2];
				
				systemState.force[i].x += pairForces[npair].x;
				systemState.force[i].y += pairForces[npair].y;
				systemState.force[i].z += pairForces[npair].z;
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
		if(step % 100 == 0 && step != 0)
			cout << "Progress: \t" << step << "\t/ " << So + Sd << endl;
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

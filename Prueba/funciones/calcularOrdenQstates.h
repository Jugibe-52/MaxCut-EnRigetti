# include <iostream>
# include "../estructuras/operation.h"
using namespace std;

void valorQStates(operation* puertas, int numeroQStates, int numeroPuertas, int** qStatesOrd);
void ordArray(int** miArray, int tamagnoArray);


int* calcularOrdenQstates(operation* puertas, int numeroQStates, int numeroPuertas){
	int** qStatesOrd;
	qStatesOrd = new int*[numeroQStates];
	int* aRetornar;
	aRetornar = new int[numeroQStates];
	valorQStates(puertas, numeroQStates, numeroPuertas, qStatesOrd);
	ordArray(qStatesOrd, numeroQStates);
	for(int i=0; i < numeroQStates; ++i){
		aRetornar[i] = qStatesOrd[i][0];
	}
	return aRetornar;
}

void valorQStates(operation* puertas, int numeroQStates, int numeroPuertas, int** qStatesOrd){
	for(int i=0; i<numeroQStates; ++i){
		qStatesOrd[i] = new int[2];
		qStatesOrd[i][0] = i;
		qStatesOrd[i][1] = 0;
	}
	for(int i=0; i < numeroPuertas; ++i){
		++qStatesOrd[puertas[i].q2][1];
		++qStatesOrd[puertas[i].q1][1];
	}
}

void ordArray(int** miArray, int tamagnoArray){
	int* temporal;
	for (int i = 0;i < tamagnoArray; i++)
	{
		for (int j = 0; j< tamagnoArray-1; j++)
		{
			if (miArray[j][1] > miArray[j+1][1])
			{
				temporal = miArray[j]; 
				miArray[j] = miArray[j+1]; 
				miArray[j+1] = temporal;
			}
		}
	}
}



# include <iostream>
using namespace std;



int* calcularFuerza(int** distHop, int dim){
	int* fuerzaQubits;
	int contador;
	fuerzaQubits = new int[dim];
	for(int i = 0; i < dim; i++){
		contador = 0;
		for(int j=0; j < dim; j++)
		{
			if(distHop[i][j] == 1 || distHop[i][j] == 2){
				contador ++;
			} 
		}
		fuerzaQubits[i] = contador;
	}
	return fuerzaQubits;
}


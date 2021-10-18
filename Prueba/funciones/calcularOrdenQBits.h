# include <iostream>
using namespace std;

void ordArrayInvert(int** miArray, int tamagnoArray);

int** calcularOrdenQBits(int* fuerzas, int numeroQBits){
	int** qBitsOrd;
	qBitsOrd = new int*[numeroQBits];
	for(int i = 0; i < numeroQBits; ++i){
		qBitsOrd[i] = new int[2];
		qBitsOrd[i][0] = i;
		qBitsOrd[i][1] = fuerzas[i];	
	}
	ordArrayInvert(qBitsOrd, numeroQBits);
	return qBitsOrd;
}

void ordArrayInvert(int** miArray, int tamagnoArray){
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

#include <iostream>
using namespace std;

// Struct representing an operation 
struct operation {
	string type; //p-s, swap or mix
	int orden = -1;
	int pos = -1;
	int posCrh = -1;
	int q1 = -1;
	int q2 = -1;
	int qb1 = -1;
	int qb2 = -1;
	int init = -1;
	int time = -1;
	int tail = -1;
	int head = -1;
	//predecesores y sucesores
	int pre1 = -1;
	int pre2 = -1;
	int suc1 = -1;
	int suc2 = -1;
};




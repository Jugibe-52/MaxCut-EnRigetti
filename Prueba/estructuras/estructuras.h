// Puertas a ejecutar 
struct puertasEjecutar {
	vector<int> a;
	vector<int> b;
	int numero;
};

// Struct representing an edge 
struct edge {
	int pos1 = -1;
	int pos2 = -1;
};

struct Individuo {
	int makespan = -1;
	int* vectorSecuenciaOperaciones;
	int* vectorPosicionOperaciones;
	int* vectorSecuenciaEdges;
	int* iniciales;
	int* qbitsBase;
	int* qstatesBase;
	int operations = -1;
	int qstatesBaseSuma = -1;
};

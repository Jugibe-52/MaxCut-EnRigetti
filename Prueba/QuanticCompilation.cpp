#include <sstream>
#include <iostream>
#include <algorithm>    // std::max
#include <time.h>
#include <math.h>
#include <vector> 
#include <windows.h>
#include <tchar.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "constantes/chipConstantes.h"
#include "constantes/constantes.h"
#include "funciones/calcularOrdenQBits.h"
#include "funciones/fuerzaConexion.h"
#include "estructuras/estructuras.h"
#include "clases/ListaQState.h"
#include "clases/ListaQBits.h"

using namespace std;

// Nombre de los parametros iniciales
const string paramIniciales[7] = {{"fichero"},{"numPruebas"},{"tiempoLimite"},{"iterLimite"},{"randomSeed"},{"poblacion"},{"iteraciones"}};

// Declaracion de los parametros iniciales
const char* rutaFichero;
int NumeroPruebas;
time_t TiempoLimite;
int IteracionesLimite;
int TAMANOPOBLACIONGENETICO;
int ITERATIONS;

// variable con las puertas que se deben ejecutar
struct puertasEjecutar puertas;

// Numero de qbits y puertas fisicas disponibles
int N, n;

int* vectorInicios;
int primero;
int ISLANDS = NUMISLANDS;      // La poblacion se divide en ''ISLANDS'' subpoblaciones, y solo los elementos de una subpoblacion se pueden cruzar entre ellos


Individuo* SolucionesIniciales;      //Para scattersearch
Individuo* Poblacion;      //Tanto para scattersearch como para genetico
Individuo* PoblacionAux;   //Para el genetico
int* parejasAleatorias;

//Datos del problema
int numTrabajos;
int numTareas;
int numMaquinas;
int** distHop;  // Distancias entre qbits
int** timePS;   // Tiempo ejecucion puertas PS
int** nextHop;  
vector<int>* sharePS;
vector<int>* allAdy;
vector<int>** edgeAdy; // puertas que son adyacentes
vector<int>** bestHop;
int** ady;
int* DuracionTareas;
int* ReleaseDates;
int* DueDates;
int* Pesos;
int** SetupTimes;                  //================SETUPS!===============================

int* ultimaIteracionQueCambio;

int numVecinos;
int cotaInferior;

//Utilizadas por la busqueda tabu
int longitudListaTabu;                          //Longitud de la lista tabu
int **cycleWitness;                             //Matriz de prevencion de ciclos
int **arcosInvertidos;                          //Matriz, cada fila es un vecino, para cada vecino hay x bloques de dos numeros, cada dos numeros indican un arco (tarea tarea) que se invierte en dicho vecino
Individuo mejorIndividuoBusquedaTabu;
Individuo indActualBusquedaTabu;
Individuo individuoaux;
Individuo individuoaux2;

int **vecinos;
int *estimaciones;
int guardados;

int mejorMakespanEncontrado;
int MejorMakespanTotalTodasPruebas = MUYPOSITIVO;
Individuo mejorIndividuoEncontrado;
Individuo MejorIndividuoTotal;
int numeroiteraciones;


int iteracionenqueseconsigueelmejor;
int iteracionenqueseconsigueelmejorOhaydiversificacion;

//Utilizadas por el cruce JOX
int* aleatorios;
int* trabajosElegidos;
int* posElegidos;
int itActual;
int itFinal;

// time_t TiempoLimite;
// int IteracionesLimite;
time_t ini;
time_t fin;

int** MatrizDistanciasSolucionesIniciales;
int** MatrizDistancias;
int distanciasBackup[TAMANOPOBLACION];
int distanciasAsiBackup[TAMANOPOBLACION];
Individuo* Trayectoria;
int MAXTRAYECTORIA;
int longitudTrayectoria;
int distanciasDeVecinosAInd2[MAXVECINOS];

int* makespanesPoblacion;
int* iteracionPoblacion;

int seleccionaux1;
int seleccionaux2;

//Para medir eficacia de recontruir cromosoma
float rcMejora = 0, rcEmpeora = 0;


// Cost of gates
#define TIME_SWAP 2
#define TIME_MIX 1


//Informacion de problema
operation* puertasDelProblema; edge* aristasDelChip;

//Qstate counters
int* qstatesTotal;
int* qstatesCount;
int* countCurrent;
int* countTotal;
int len;

//contar vecinos

float vec1T = 0;
float vec2T = 0;
float vec3T = 0;
float vec1 = 0;
float vec2 = 0;
float vec3 = 0;
float vec1m = 0;
float vec2m = 0;
float vec3m = 0;
float vec1e = 0;
float vec2e = 0;
float vec3e = 0;
float bloqueos = 0;


//Para la segunda iteracion
vector<vector<operation>> mejoresInicios;
bool guardar;

/*METODOS*/
void IniciarArgumentos(const char** paramProb);
struct puertasEjecutar LeerPuertas(const char * nombre);
void InicioChip();
void IniciarQStates();
void create_gnome(vector<int> q1, vector<int> q2, int len, operation* genOp, int iterations, int* qTotal, int* countCurrent);
int Aleatorio(int minimo, int maximo);
void IniciarIndividuo(Individuo &ind);
void InicializaEstructuras(void);
void CrearVectoresIniciales(Individuo &ind);
void PlanificacionSemiactiva(Individuo &ind, int itIni, int itFin);
void PlanificacionSemiactiva(Individuo &ind, vector<operation> *operations, int itIni, int itFin);
void CopiarIndividuo(Individuo &ind1, Individuo &ind2, int nivel);
void CrearPoblacionInicialGenetico(void);
void Cruce(int ind1, int ind2);
void Mutacion(Individuo &ind, int c, int t);
bool Iguales(Individuo &ind1, Individuo &ind2);
void BorraCabezasYColas(vector<operation>* operations);
void insertSWAP(edge ed, operation *op, vector<operation>* operations, int *qbits, int len, int *qstates);
void insertPS(operation *op, vector<operation>* operations, int *qbits, int *qstates, int len, int* qstatesTotal, int* qstatesCount);

//Declaraciones
void addSwap(int i, int si, edge &ed, vector<operation>* operations, int *qbits, int *qstates);
void LocalizarQstates(int * qbits, int * qstates, vector<operation> *operations, operation *op, int *t1, int *i1, int *t2, int *i2);
void addMix(int * qstatesCount, int qmix, int * qstatesTotal, int t, vector<operation>* operations, int i, int* qstates);
void GeneraSchedule(std::vector<operation> *operations, Individuo & ind, int it, int it2);
void PrintSchedule(std::vector<operation> &operations1);
void PrintSchedule2(std::vector<operation> &operations1);

// Generar Vector Aleatorio TFM
int  existeElem (int* v, int longi, int num);
void generarVector (int* v);
void escribirVec (int* v, int longi);
void restarUno (int* v, int longi);

void insertarQBitLista(ListaQBits* listaQBitsInic, ListaQBits* listaQBitsFinal, QState* sigQState);

//FUNCION PRINCIPAL
int main(int argc,const char** argv) {
	// Codigo comentado hasta que lo comprenda
	// SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
	if (argc < 8) {
		cout << endl << "Uso: QuanticCompilation fichero numpruebas tiempolimite iteracioneslimite randomseed population iterations" << endl;
		return 0;
	}
	
	// Imprimimos los parametros iniciales
	for (int i = 0 ; i < argc - 1 ; ++i)
		cout << paramIniciales[i] << ": "<< argv[i + 1] << endl;
		
	cout << endl << endl;	
	
	// Iniciamos los argumentos que entran por consola.
	IniciarArgumentos(argv);
	
	// Leemos y guardamos las puertas a ejecutar
	puertas = LeerPuertas(rutaFichero);
		
	// Iniciamos chip (variables distHop, )
	InicioChip();
	
	// Iniciamos QStates
	IniciarQStates();

	//Iniciza estructuras
	InicializaEstructuras();
	


	int **parejasHechas = new int*[TAMANOPOBLACION];
	for (int i = 0; i < TAMANOPOBLACION; i++) parejasHechas[i] = new int[TAMANOPOBLACION];
	
	if (METODO == 1) ultimaIteracionQueCambio = new int[TAMANOPOBLACIONGENETICO];
	else ultimaIteracionQueCambio = new int[TAMANOPOBLACION];
	
	int* makespanPruebaNumero = new int[NumeroPruebas];
	int* firstTimePruebaNumero = new int[NumeroPruebas];
	int* firstIterPruebaNumero = new int[NumeroPruebas];
	
	// Logs finales
	int NumeroMejorMakespanesTotalesTodasPruebas = 0;
	int PeorMakespanTotalTodasPruebas = MUYNEGATIVO;
	int NumeroPeoresMakespanesTotalesTodasPruebas = 0;
	double MediaMakespanTotalTodasPruebas = 0;
	double DesviacionEstandarTotalTodasPruebas = 0;
	double TiempoTotalTodasPruebas = 0;

	int totalIteracionesTodasPruebas = 0;
	int totalDiversificaciones = 0;
	
	

	// (TFM) Inicializacion de qstates, guardar, copiar mas adelante.
	
	// Iniciar vectorInicio
    int longi=N, num;
	int numeroQStates = puertas.numero;
	int* qStatesOrdInvert = calcularOrdenQstates(puertasDelProblema, numeroQStates, N);
	int* fuerzaQBits = calcularFuerza(distHop, N);
	int** qBitsOrdInvert = calcularOrdenQBits(fuerzaQBits, N);
	
	ListaQBits* listaQBitsInic = new ListaQBits();
	listaQBitsInic->iniciarListaConFuer(qBitsOrdInvert, N);
	
	listaQBitsInic->imprimirFuerza();
	
	ListaQState* listaQStatesInic = new ListaQState();
	listaQStatesInic->iniciarLista(qStatesOrdInvert, numeroQStates);
	
	ListaQBits* listaQBitsFinal = new ListaQBits();

	listaQStatesInic->imprimir();
	
	QState* sigQState;
	QBit* sigQBit;
	sigQState = listaQStatesInic->extraerPrimero();
	
	cout << "Comenzamos la inicializacion de los QStates:" << endl << endl;
	insertarQBitLista(listaQBitsInic, listaQBitsFinal, sigQState);

	
	while(listaQStatesInic->getRaiz() != NULL)
	{
		sigQState = listaQStatesInic->extraerPrimero();
		cout << "INICIALIZACION DEL QState: " << sigQState->getQState() + 1<< endl << endl;
		sigQBit = listaQBitsInic->extraerHeredero(listaQBitsInic, listaQBitsFinal, sigQState,  puertasDelProblema, N, distHop);
		if(sigQBit != NULL)
		{
			sigQBit = listaQBitsInic->extraerQBit(sigQBit);
			cout << "Insertamos el QState " << sigQBit->qStateCont->getQState() + 1<< " en el QBit " << sigQBit->getNumero() + 1 << endl;
			listaQBitsFinal->insertarQBit(sigQBit);
		}
		else
		{
			cout << "No se han insertado Vecinos logicos para el QState " << sigQState->getQState() + 1 << endl;
			insertarQBitLista(listaQBitsInic, listaQBitsFinal, sigQState);	
		}
		cout << endl;
	}

	
	int qStatesAcomp = numeroQStates;
	while(listaQBitsInic->getRaiz() != NULL)
	{
		QBit* qBitAIns = listaQBitsInic->extraerPrimero();
		QState* extrQState=new QState(qStatesAcomp);
		qBitAIns->qStateCont = extrQState;
		listaQBitsFinal->insertarQBit(qBitAIns);
		qStatesAcomp++;
		
	}
	
	listaQBitsFinal->imprimir();
	vectorInicios = new int[N];

	listaQBitsFinal->iniciarLista(vectorInicios);
	for(int i = 0; i<N; ++i){
		vectorInicios[i] = listaQBitsFinal->devolverState(i);
	}
	


	for (int numpru = 0; numpru < NumeroPruebas; numpru++) {

		ini = time(NULL);
		fin = time(NULL);
		primero = PRIMEROS;  //Numero de individuos que se generan de forma heuritica al crear la poblacion
		ISLANDS = NUMISLANDS;  // La poblacion se divide en ''ISLANDS'' subpoblaciones, y solo los elementos de una subpoblacion se pueden cruzar entre ellos
		guardar = false;

		int numdiversificaciones = 0;

		for (int i = 0; i < TAMANOPOBLACION - 1; i++) for (int j = i + 1; j < TAMANOPOBLACION; j++) parejasHechas[i][j] = -1;


		fin = time(NULL);


		if (METODO == 1) {        //ALGORITMO GENETICO GENERACIONAL

			itActual = 0;
			if (ETAPAS == 1)      //Si esta a 1 realiza la busqueda por etapas, tratando solo una iteracion de cada vez
				itFinal = 1;
			else
				itFinal = ITERATIONS;
			mejorMakespanEncontrado = MUYPOSITIVO;
			
			CrearPoblacionInicialGenetico();


			for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) 
			{
				if (Poblacion[i].makespan < mejorMakespanEncontrado) {
					mejorMakespanEncontrado = Poblacion[i].makespan;
					CopiarIndividuo(Poblacion[i], mejorIndividuoEncontrado, 1);
					iteracionenqueseconsigueelmejor = numeroiteraciones;
					iteracionenqueseconsigueelmejorOhaydiversificacion = numeroiteraciones;
				}
				ultimaIteracionQueCambio[i] = numeroiteraciones;
			}

			fin = time(NULL);
			numeroiteraciones = 0;
			iteracionenqueseconsigueelmejor = 0;
			iteracionenqueseconsigueelmejorOhaydiversificacion = 0;
			int probabilidadmutacion = PROBMUTACION;

			mejoresInicios = vector<vector<operation>>();

			//Bucle principal
			while (((TiempoLimite == 0 || (TiempoLimite != 0 && (fin - ini < TiempoLimite))) && numeroiteraciones - iteracionenqueseconsigueelmejor < MAXGENERACIONESSINMEJORAGENETICO) || itFinal != ITERATIONS || ETAPASF == 1 && itActual != 0) {

				if (mejorMakespanEncontrado <= cotaInferior)
					break;

				//Genera la poblacion inicial cuando se cambia de fase
				if (numeroiteraciones - iteracionenqueseconsigueelmejor >= MAXGENERACIONESSINMEJORAGENETICO) {

					guardados = 0;
					guardar = true;
					for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
						PlanificacionSemiactiva(Poblacion[i], itActual, itFinal);
					}

					PlanificacionSemiactiva(mejorIndividuoEncontrado, itActual, itFinal);

					//Y actualizar las variables necesarias
					for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
						

						if (Poblacion[i].makespan < mejorMakespanEncontrado) {
							mejorMakespanEncontrado = Poblacion[i].makespan;
							CopiarIndividuo(Poblacion[i], mejorIndividuoEncontrado, 1);
						}
						ultimaIteracionQueCambio[i] = numeroiteraciones;
					}

					firstTimePruebaNumero[numpru] = fin - ini;
					firstIterPruebaNumero[numpru] = numeroiteraciones;

					itActual++;
					itFinal++;
					if (MUTACIONADAPTATIVA == 1) probabilidadmutacion = 0;
					primero = PRIMEROS;


					//Actualizar iteracion
					if (itActual == ITERATIONS) {
						itActual = 0;
						itFinal = ITERATIONS;
					}

					guardar = false;

					
					//Crear nueva poblacion
					CrearPoblacionInicialGenetico();

					//Y actualizar las variables necesarias
					for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {

						if (Poblacion[i].makespan < mejorMakespanEncontrado) {
							mejorMakespanEncontrado = Poblacion[i].makespan;
							CopiarIndividuo(Poblacion[i], mejorIndividuoEncontrado, 1);
							iteracionenqueseconsigueelmejor = numeroiteraciones;
							iteracionenqueseconsigueelmejorOhaydiversificacion = numeroiteraciones;
						}
						ultimaIteracionQueCambio[i] = numeroiteraciones;
					}

					iteracionenqueseconsigueelmejorOhaydiversificacion = numeroiteraciones;
				}


				//AGRUPAR LA POBLACION ALEATORIAMENTE POR PAREJAS
				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) parejasAleatorias[i] = i;
				for (int i = 0; i < ISLANDS; i++) 
				{
					int inicioisla = (TAMANOPOBLACIONGENETICO / ISLANDS)*i;
					int finisla = inicioisla + (TAMANOPOBLACIONGENETICO / ISLANDS);
					for (int j = inicioisla; j < finisla; j++) 
					{
						int nuevaPosicionAleatoria = Aleatorio(inicioisla, finisla - 1);
						int posicionbackup = parejasAleatorias[j];
						parejasAleatorias[j] = parejasAleatorias[nuevaPosicionAleatoria];
						parejasAleatorias[nuevaPosicionAleatoria] = posicionbackup;
					}
				}

				//APLICAR EL OPERADOR DE CRUCE (pilla Poblacion[ind1] y Poblacion[ind2] y los cruza en PoblacionAux[ind1] y PoblacionAux[ind2])
				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i = i + 2) {
					for (int j = 0; j < numTareas; j++) {
						CopiarIndividuo(Poblacion[parejasAleatorias[i]], PoblacionAux[parejasAleatorias[i]], 1);
						CopiarIndividuo(Poblacion[parejasAleatorias[i + 1]], PoblacionAux[parejasAleatorias[i + 1]], 1);
					}
					int aleat = Aleatorio(1, 100);
					if (aleat <= PROBCRUCE) Cruce(parejasAleatorias[i], parejasAleatorias[i + 1]);
				}

				//APLICAR EL OPERADOR DE MUTACION
				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
					int aleat = Aleatorio(1, 100);
					if (aleat <= probabilidadmutacion) Mutacion(PoblacionAux[i], TIPOMUTACION, 0);
				}

				//APLICAR LA BUSQUEDA TABU A TODOS
				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
					PlanificacionSemiactiva(PoblacionAux[i], itActual, itFinal);
				}

				if (DEBUG >= 2) {
					cout << endl << "POBLACION AL PRINCIPIO DE LA GENERACION " << numeroiteraciones << ":" << endl;
					for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) cout << Poblacion[i].makespan << " ";
					cout << endl;
					if (DEBUG >= 3) {
						cout << endl << "PAREJASALEATORIAS: ";
						for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) cout << parejasAleatorias[i] << " ";
						cout << endl << "POBLACION DESPUES DE CRUZAR Y APLICAR BUSQUEDAS TABU: " << endl;
						for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) cout << PoblacionAux[i].makespan << endl;
					}
					char c;
					cin >> c;
				}

				//ELEGIR LOS DOS MEJORES DE CADA PAREJA DE PADRES E HIJOS
				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i = i + 2) 
				{
					int pos1 = parejasAleatorias[i];
					int pos2 = parejasAleatorias[i + 1];
					int make1 = Poblacion[pos1].makespan;
					int make2 = Poblacion[pos2].makespan;
					int make3 = PoblacionAux[pos1].makespan;
					int make4 = PoblacionAux[pos2].makespan;

					//cout << endl << make1 << " " << make2 << " " << make3 << " " << make4;

					if (SELECCION == 1) {  //Elegir los dos mejores
						if ((make3 <= make1) && (make3 <= make2) && (make4 <= make1) && (make4 <= make2)) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
						}
						else if ((make1 <= make2) && (make1 <= make4) && (make3 <= make2) && (make3 <= make4)) CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos2], 1);
						else if ((make1 <= make2) && (make1 <= make3) && (make4 <= make2) && (make4 <= make3)) CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
						else if ((make2 <= make1) && (make2 <= make4) && (make3 <= make1) && (make3 <= make4)) CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
						else if ((make2 <= make1) && (make2 <= make3) && (make4 <= make1) && (make4 <= make3)) CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos1], 1);
					}

					if (SELECCION == 2) {  //Elegir el mejor y el peor
						int mejorposs = 1;
						int mejorr = make1;
						if (make2 < mejorr) { mejorposs = 2; mejorr = make2; }
						if (make3 < mejorr) { mejorposs = 3; mejorr = make3; }
						if (make4 < mejorr) { mejorposs = 4; }
						int peorposs = 4;
						int peorr = make4;
						if (make3 > peorr) { peorposs = 3; peorr = make3; }
						if (make2 > peorr) { peorposs = 2; peorr = make2; }
						if (make1 > peorr) { peorposs = 1; }
						if (((mejorposs == 3) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 3))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
						}
						else if (((mejorposs == 1) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 1))) CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos2], 1);
						else if (((mejorposs == 1) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 1))) CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
						else if (((mejorposs == 2) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 2))) CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
						else if (((mejorposs == 2) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 2))) CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos1], 1);
					}

					if (SELECCION == 3) {  //Elegir los dos mejores pero tales que tengan diferente makespan
						int mejorposs = 4;
						int mejorr = make4;
						if (make3 < mejorr) { mejorposs = 3; mejorr = make3; }
						if (make2 < mejorr) { mejorposs = 2; mejorr = make2; }
						if (make1 < mejorr) { mejorposs = 1; mejorr = make1; }
						int peorposs;
						int peorr;
						seleccionaux1++;
						if ((make1 == make2) && (make2 == make3) && (make3 == make4)) {
							peorposs = 3;
						}
						else {
							if (mejorposs == 1) {
								peorposs = 4;
								peorr = make4;
								if (((make3 < peorr) && (make3 != mejorr)) || (peorr == mejorr)) { peorposs = 3; peorr = make3; }
								if (((make2 < peorr) && (make2 != mejorr)) || (peorr == mejorr)) { peorposs = 2; }
							}
							else if (mejorposs == 2) {
								peorposs = 4;
								peorr = make4;
								if (((make3 < peorr) && (make3 != mejorr)) || (peorr == mejorr)) { peorposs = 3; peorr = make3; }
								if (((make1 < peorr) && (make1 != mejorr)) || (peorr == mejorr)) { peorposs = 1; }
							}
							else if (mejorposs == 3) {
								peorposs = 4;
								peorr = make4;
								if (((make2 < peorr) && (make2 != mejorr)) || (peorr == mejorr)) { peorposs = 2; peorr = make2; }
								if (((make1 < peorr) && (make1 != mejorr)) || (peorr == mejorr)) { peorposs = 1; }
							}
							else if (mejorposs == 4) {
								peorposs = 3;
								peorr = make3;
								if (((make2 < peorr) && (make2 != mejorr)) || (peorr == mejorr)) { peorposs = 2; peorr = make2; }
								if (((make1 < peorr) && (make1 != mejorr)) || (peorr == mejorr)) { peorposs = 1; }
							}
						}

						if (((mejorposs == 3) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 3))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
							//	          cout << " " << "3 y 4";
							//			  if((make1<make3)||(make1<make4)||(make2<make3)||(make2<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 1) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 1))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos2], 1);
							//			  cout << " " << "1 y 3";
							//			  if((make2<make1)||(make2<make3)||(make4<make1)||(make4<make3)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 1) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 1))) {
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
							//		      cout << " " << "1 y 4";
							//			  if((make2<make1)||(make2<make4)||(make3<make1)||(make3<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 2) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 2))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							//			  cout << " " << "2 y 3";
							//			  if((make1<make2)||(make1<make3)||(make4<make2)||(make4<make3)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 2) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 2))) {
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos1], 1);
							//			  cout << " " << "2 y 4";
							//			  if((make1<make2)||(make1<make4)||(make3<make2)||(make3<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						//			else if (((mejorposs==1)&&(peorposs==2))||((mejorposs==2)&&(peorposs==1))) {
						//			  cout << " " << "1 y 2";
						//			  if((make3<make1)||(make3<make2)||(make4<make1)||(make4<make2)) { seleccionaux2++; cout << " EJEM"; }
						//			}

					}

					if (SELECCION == 4) {  //Elegir los dos mejores pero tales que sean diferentes
		  //cout << endl << endl << make1 << " " << make2 << " " << make3 << " " << make4 << endl;
		  //for(int ghj=0; ghj<numTareas; ghj++) cout << Poblacion[pos1].vectorSecuencia[ghj] << " ";
		  //cout << endl;
		  //for(int ghj=0; ghj<numTareas; ghj++) cout << Poblacion[pos2].vectorSecuencia[ghj] << " ";
		  //cout << endl;
		  //for(int ghj=0; ghj<numTareas; ghj++) cout << PoblacionAux[pos1].vectorSecuencia[ghj] << " ";
		  //cout << endl;
		  //for(int ghj=0; ghj<numTareas; ghj++) cout << PoblacionAux[pos2].vectorSecuencia[ghj] << " ";
		  //cout << endl;
						int mejorposs = 4;
						int mejorr = make4;
						if (make3 < mejorr) { mejorposs = 3; mejorr = make3; }
						if (make2 < mejorr) { mejorposs = 2; mejorr = make2; }
						if (make1 < mejorr) { mejorposs = 1; mejorr = make1; }
						int peorposs;
						int peorr;
						seleccionaux1++;
						if ((Iguales(Poblacion[pos1], Poblacion[pos2]) == true) && (Iguales(Poblacion[pos2], PoblacionAux[pos1]) == true) && (Iguales(PoblacionAux[pos1], PoblacionAux[pos2]) == true)) {
							peorposs = 3;
						}
						else {
							if (mejorposs == 1) {
								peorposs = 4;
								peorr = make4;
								if (((make3 < peorr) && (Iguales(Poblacion[pos1], PoblacionAux[pos1]) == false)) || (Iguales(Poblacion[pos1], PoblacionAux[pos2]) == true)) { peorposs = 3; peorr = make3; }
								if (peorposs == 4) {
									if (((make2 < peorr) && (Iguales(Poblacion[pos2], Poblacion[pos1]) == false)) || (Iguales(PoblacionAux[pos2], Poblacion[pos1]) == true)) { peorposs = 2; }
								}
								else {
									if (((make2 < peorr) && (Iguales(Poblacion[pos2], Poblacion[pos1]) == false)) || (Iguales(PoblacionAux[pos1], Poblacion[pos1]) == true)) { peorposs = 2; }
								}
							}
							else if (mejorposs == 2) {
								peorposs = 4;
								peorr = make4;
								if (((make3 < peorr) && (Iguales(PoblacionAux[pos1], Poblacion[pos2]) == false)) || (Iguales(Poblacion[pos2], PoblacionAux[pos2]) == true)) { peorposs = 3; peorr = make3; }
								if (peorposs == 4) {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], Poblacion[pos2]) == false)) || (Iguales(PoblacionAux[pos2], Poblacion[pos2]) == true)) { peorposs = 1; }
								}
								else {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], Poblacion[pos2]) == false)) || (Iguales(PoblacionAux[pos1], Poblacion[pos2]) == true)) { peorposs = 1; }
								}
							}
							else if (mejorposs == 3) {
								peorposs = 4;
								peorr = make4;
								if (((make2 < peorr) && (Iguales(Poblacion[pos2], PoblacionAux[pos1]) == false)) || (Iguales(PoblacionAux[pos1], PoblacionAux[pos2]) == true)) { peorposs = 2; peorr = make2; }
								if (peorposs == 4) {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], PoblacionAux[pos1]) == false)) || (Iguales(PoblacionAux[pos2], PoblacionAux[pos1]) == true)) { peorposs = 1; }
								}
								else {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], PoblacionAux[pos1]) == false)) || (Iguales(Poblacion[pos2], PoblacionAux[pos1]) == true)) { peorposs = 1; }
								}
							}
							else if (mejorposs == 4) {
								peorposs = 3;
								peorr = make3;
								if (((make2 < peorr) && (Iguales(Poblacion[pos2], PoblacionAux[pos2]) == false)) || (Iguales(PoblacionAux[pos1], PoblacionAux[pos2]) == true)) { peorposs = 2; peorr = make2; }
								if (peorposs == 3) {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], PoblacionAux[pos2]) == false)) || (Iguales(PoblacionAux[pos1], PoblacionAux[pos2]) == true)) { peorposs = 1; }
								}
								else {
									if (((make1 < peorr) && (Iguales(Poblacion[pos1], PoblacionAux[pos2]) == false)) || (Iguales(Poblacion[pos2], PoblacionAux[pos2]) == true)) { peorposs = 1; }
								}
							}
						}

						if (((mejorposs == 3) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 3))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
							//            cout << " " << "3 y 4";
							//			  if((make1<make3)||(make1<make4)||(make2<make3)||(make2<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 1) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 1))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos2], 1);
							//            cout << " " << "1 y 3";
							//			  if((make2<make1)||(make2<make3)||(make4<make1)||(make4<make3)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 1) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 1))) {
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos2], 1);
							//            cout << " " << "1 y 4";
							//			  if((make2<make1)||(make2<make4)||(make3<make1)||(make3<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 2) && (peorposs == 3)) || ((mejorposs == 3) && (peorposs == 2))) {
							CopiarIndividuo(PoblacionAux[pos1], Poblacion[pos1], 1);
							//            cout << " " << "2 y 3";
							//			  if((make1<make2)||(make1<make3)||(make4<make2)||(make4<make3)) { seleccionaux2++; cout << " EJEM"; }
						}
						else if (((mejorposs == 2) && (peorposs == 4)) || ((mejorposs == 4) && (peorposs == 2))) {
							CopiarIndividuo(PoblacionAux[pos2], Poblacion[pos1], 1);
							//            cout << " " << "2 y 4";
							//			  if((make1<make2)||(make1<make4)||(make3<make2)||(make3<make4)) { seleccionaux2++; cout << " EJEM"; }
						}
						//          else if (((mejorposs==1)&&(peorposs==2))||((mejorposs==2)&&(peorposs==1))) {
						//            cout << " " << "1 y 2";
						//			  if((make3<make1)||(make3<make2)||(make4<make1)||(make4<make2)) { seleccionaux2++; cout << " EJEM"; }
						//          }
						//char c;
						//cin >> c;
					}


				}

				//ACTUALIZAR iteracionenqueseconsigueelmejor Y EL CONTADOR DE GENERACIONES
				numeroiteraciones++;

				for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {

					if (Poblacion[i].makespan < mejorMakespanEncontrado) {
						mejorMakespanEncontrado = Poblacion[i].makespan;
						CopiarIndividuo(Poblacion[i], mejorIndividuoEncontrado, 1);
						iteracionenqueseconsigueelmejor = numeroiteraciones;
						iteracionenqueseconsigueelmejorOhaydiversificacion = numeroiteraciones;
					}
					ultimaIteracionQueCambio[i] = numeroiteraciones;
				}

				//APLICAR ALGO DE DIVERSIFICACION
				//Realiza algunas mutacones en los individuos con makespan repetido
				if ((ELIMINARIGUALES > 0) && (numeroiteraciones > 0) && (numeroiteraciones - iteracionenqueseconsigueelmejor - 2) >= ELIMINARIGUALES && (((numeroiteraciones - iteracionenqueseconsigueelmejor - 2) % ELIMINARIGUALES) == 0)) {

					for (int i = 0; i < ISLANDS; i++) {
						int inicioisla = (TAMANOPOBLACIONGENETICO / ISLANDS)*i;
						int finisla = inicioisla + (TAMANOPOBLACIONGENETICO / ISLANDS);

						for (int k = inicioisla; k < finisla; k++) {
							if (iteracionPoblacion[Poblacion[k].makespan] == (numeroiteraciones*ISLANDS) + i) {
								makespanesPoblacion[Poblacion[k].makespan]++;
							}
							else {
								iteracionPoblacion[Poblacion[k].makespan] = (numeroiteraciones*ISLANDS) + i;
								makespanesPoblacion[Poblacion[k].makespan] = 1;
							}
						}
						for (int k = inicioisla; k < finisla; k++) {
							if ((makespanesPoblacion[Poblacion[k].makespan] > 1) || (Poblacion[k].makespan >= (mejorMakespanEncontrado + 1) && makespanesPoblacion[Poblacion[k].makespan] > 1)) {
								makespanesPoblacion[Poblacion[k].makespan]--;
								if (Aleatorio(1, 100) <= ELIMINARIGUALESPOR) Mutacion(Poblacion[k], Aleatorio(1, n / 2) , 1);
								
								PlanificacionSemiactiva(Poblacion[k], itActual, itFinal);
							}
						}
					}
				}


				//Y tratar la mutacion adaptativa si es necesario
				if (MUTACIONADAPTATIVA == 1) {
					if (iteracionenqueseconsigueelmejor == numeroiteraciones) probabilidadmutacion = 0;
					else if (probabilidadmutacion < PROBMUTACION) probabilidadmutacion = probabilidadmutacion + PORCENTMUTACIONADAPTATIVA;
				}


				//Informacion de debug
				if (DEBUG >= 1) {
					int bestm = MUYPOSITIVO;
					double avgm = 0;
					for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
						if (Poblacion[i].makespan < bestm) bestm = Poblacion[i].makespan;
						avgm = avgm + Poblacion[i].makespan;
					}
					avgm = avgm / TAMANOPOBLACIONGENETICO;
					
					//if (printIteration = true)
					//cout << endl << numeroiteraciones << " " << bestm << " " << avgm;
				}

				fin = time(NULL);
			}
			vector<vector<operation>>().swap(mejoresInicios);
		}


		//Actualiza el mejor makesman e individuo encontrado
		makespanPruebaNumero[numpru] = mejorMakespanEncontrado;
		if (mejorMakespanEncontrado < MejorMakespanTotalTodasPruebas) {
			MejorMakespanTotalTodasPruebas = mejorMakespanEncontrado;
			CopiarIndividuo(mejorIndividuoEncontrado, MejorIndividuoTotal, 1);
			NumeroMejorMakespanesTotalesTodasPruebas = 1;
		}
		else if (mejorMakespanEncontrado == MejorMakespanTotalTodasPruebas) NumeroMejorMakespanesTotalesTodasPruebas++;
		if (mejorMakespanEncontrado > PeorMakespanTotalTodasPruebas) {
			PeorMakespanTotalTodasPruebas = mejorMakespanEncontrado;
			NumeroPeoresMakespanesTotalesTodasPruebas = 1;
		}
		else if (mejorMakespanEncontrado == PeorMakespanTotalTodasPruebas) NumeroPeoresMakespanesTotalesTodasPruebas++;
		MediaMakespanTotalTodasPruebas += mejorMakespanEncontrado;
		TiempoTotalTodasPruebas += (fin - ini);

		totalIteracionesTodasPruebas += numeroiteraciones;
		totalDiversificaciones += numdiversificaciones;

	}


	//Calcula estadisticas de las pruebas
	int mediana = 0; int nummediana = 0; int firsTime = 0; int firstIter = 0;
	MediaMakespanTotalTodasPruebas = MediaMakespanTotalTodasPruebas / NumeroPruebas;
	if (NumeroPruebas != 1) {
		double sumatorio = 0;
		for (int i = 0; i < NumeroPruebas; i++) {
			sumatorio += (((double)makespanPruebaNumero[i] - MediaMakespanTotalTodasPruebas)*((double)makespanPruebaNumero[i] - MediaMakespanTotalTodasPruebas));
		}
		DesviacionEstandarTotalTodasPruebas = sqrt(sumatorio / (double)(NumeroPruebas - 1));

		for (int i = 0; i < NumeroPruebas; i++) {
			int count = 0;
			firsTime += firstTimePruebaNumero[i];
			firstIter += firstIterPruebaNumero[i];
			for (int j = 0; j < NumeroPruebas; j++) {
				if (makespanPruebaNumero[i] == makespanPruebaNumero[j])
					count++;
			}
			if (count > nummediana) {
				mediana = makespanPruebaNumero[i];
				nummediana = count;
			}
		}

	}
	else {
		firsTime = firstTimePruebaNumero[0];
		firstIter = firstIterPruebaNumero[0];
	}

	int total = vec1 + vec2 + vec3;

	//Imprime informacion
	cout << "Mejor Makespan Total Todas Pruebas: " << MejorMakespanTotalTodasPruebas << " \n" 
	<< "Numero Mejor Makespanes Totales Todas Pruebas: " << NumeroMejorMakespanesTotalesTodasPruebas << " \n" 
	<< "Media Makespan Total Todas Pruebas: " << MediaMakespanTotalTodasPruebas << " \n" 
	<< "Mediana: " << mediana << " \n" 
	<< "Ejecucione en la mediana: " << nummediana << " \n"
	<< "Peor Makespan Total Todas Pruebas: " << PeorMakespanTotalTodasPruebas << " \n" 
	<< "Numero Peores Makespanes Totales Todas Pruebas: " << NumeroPeoresMakespanesTotalesTodasPruebas << " \n" 
	<< "Desviacion Estandar Total Todas Pruebas: " << DesviacionEstandarTotalTodasPruebas << " \n" 
	<< "Tiempo Total Todas Pruebas: " << TiempoTotalTodasPruebas << " \n";

	cout << MejorMakespanTotalTodasPruebas << "&"
		<< PeorMakespanTotalTodasPruebas << "&"
		<< MediaMakespanTotalTodasPruebas << "&"
		<< TiempoTotalTodasPruebas << "&"
		<< DesviacionEstandarTotalTodasPruebas << "&" << endl;
	vector<operation> operations = vector<operation>();
	vector<int> posOperation = vector<int>();

	//MaximoGradiente(MejorIndividuoTotal, 0, &operations);
	PlanificacionSemiactiva(MejorIndividuoTotal, &operations, 0, ITERATIONS);

	BorraCabezasYColas(&operations);
	PrintSchedule(operations);	cout << "\n";
	PrintSchedule2(operations);	cout << "\n";
	

	return 1;

}

void insertarQBitLista(ListaQBits* listaQBitsInic, ListaQBits* listaQBitsFinal, QState* sigQState){
	QBit* qBitAInsertar = listaQBitsInic->extraerPrimero();
	qBitAInsertar->qStateCont = sigQState;
	cout << "Insertamos el QState " << qBitAInsertar->qStateCont->getQState() + 1
	<< " en el QBit " << qBitAInsertar->getNumero() + 1<< endl << endl;
	listaQBitsFinal->insertarQBit(qBitAInsertar);
}

void IniciarArgumentos(const char** paramProb){
	/* Parametros de la ejecucion */
	// Ruta hacia el fichero que lee las puertas a ejecutar
	rutaFichero = paramProb[1];
	// Numero de schedulings que se generan
	NumeroPruebas = atoi(paramProb[2]);
	// Tiempo limite de la ejecucion, normalmente a cero
	TiempoLimite = atoi(paramProb[3]);
	// Numero de iteraciones realizadas para calcular la ultima generacion
	IteracionesLimite = atoi(paramProb[4]);
	// Semilla utilizada en el problema
	srand(atoi(paramProb[5]));
	// Tamano de la poblacion con cada generacion
	TAMANOPOBLACIONGENETICO = atoi(paramProb[6]);
	// Numero de veces que hay que ejecutar las puertas cu�nticas
	ITERATIONS = atoi(paramProb[7]);	
}

struct puertasEjecutar LeerPuertas(const char * nombre){
	
	struct puertasEjecutar thispuertas;
	vector<int> a = vector<int>();
	vector<int> b = vector<int>();
	fstream f(nombre);
	string aux;
	
	//Buscar datos problema
	getline(f, aux);
	if (aux == "(define (problem compiledcode)") {

		while (aux != "(:goal")
		{
			getline(f, aux);
		}
		getline(f, aux);
		getline(f, aux);
		int max = 0;
		numTrabajos = 0;
		while (aux != "   )")
		{
			numTrabajos++;
			bool found = false;
			bool fo = false;
			int j = 0;
			for (int i = 0; i < aux.length(); i++) {
				if (isdigit(aux[i])) {
					if (!found) {
						j = i;
						found = true;
					}
				}
				else if (found) {
					//almacena el numero que se guarda en a o en b
					int nl = atoi(aux.substr(j, i - j).c_str());
					found = false;
					//insert number
					if (!fo) {
						a.push_back(nl);
						fo = true;
					}
					else
						b.push_back(nl);
					if (max < nl)
						max = nl;
				}
			}
			getline(f, aux);
		}
				
		// sobrescribimos el numero de puertas
		thispuertas.numero = max;	
	}
	f.close();
	thispuertas.a = a;
	thispuertas.b = b;
	return thispuertas;
}

void InicioChip() {
	vector<int> a, b;
	int cols, rows;
	
	// Tipo chip
	if (puertas.numero > 21){
		N = 40;
		cols = 40, rows = 40;
	}
		
	else if (puertas.numero > 8){
		N = 21;
		cols = 21, rows = 21;
	}
	else if (puertas.numero > 4){
		N = 8;
		cols = 8, rows = 8;
	}	
	else{
		N = 4;
		cols = 4, rows = 4;
	}
		
	distHop = new int*[rows];
	timePS = new int*[rows];
	// nextHop = new int*[rows];
	for (int i = 0; i < rows; ++i) {
		distHop[i] = new int[cols];
		timePS[i] = new int[cols];
		// nextHop[i] = new int[cols];
		if(N == 4){
			for (int j = 0; j < cols; ++j) {
				distHop[i][j] = distHop4[i][j];
				timePS[i][j] = timePS4[i][j];
				//nextHop[i][j] = nextHop4[i][j];
			}
		} else if (N == 8) {
			for (int j = 0; j < cols; ++j) {
				distHop[i][j] = distHop8[i][j];
				timePS[i][j] = timePS8[i][j];
				// nextHop[i][j] = nextHop4[i][j];
			}
		} else if (N == 21) {
			for (int j = 0; j < cols; ++j) {
				distHop[i][j] = distHop21[i][j];
				timePS[i][j] = timePS21[i][j];
				// nextHop[i][j] = nextHop4[i][j];
			}
		} else if (N == 40) {
			for (int j = 0; j < cols; ++j) {
				distHop[i][j] = distHop40[i][j];
				timePS[i][j] = timePS40[i][j];
				// nextHop[i][j] = nextHop4[i][j];
			}
		} else if (N == 72) {
			for (int j = 0; j < cols; ++j) {
				distHop[i][j] = distHop72[i][j];
				timePS[i][j] = timePS72[i][j];
				// nextHop[i][j] = nextHop4[i][j];
			}
		}
	}
	if(N == 4){
		cols = 2; rows = 4;
	} else if (N == 8) {
		cols = 2; rows = 8;
	} else if (N == 21) {
		cols = 2; rows = 24;
	} else if (N == 40) {
		cols = 2, rows = 48;
	} else if (N == 72) {
		cols = 2, rows = 121;
	}
	ady = new int*[rows];
	for (int i = 0; i < rows; ++i) {
		ady[i] = new int[cols];
		if(N == 4){
			for (int j = 0; j < cols; ++j) {
				ady[i][j] = ady4[i][j];
			}
		} else if (N == 8) {
			for (int j = 0; j < cols; ++j) {
				ady[i][j] = ady8[i][j];
			}
		} else if (N == 21) {
			for (int j = 0; j < cols; ++j) {
				ady[i][j] = ady21[i][j];
			}
		} else if (N == 40) {
			for (int j = 0; j < cols; ++j) {
				ady[i][j] = ady40[i][j];
			}
		} else if (N == 72) {
			for (int j = 0; j < cols; ++j) {
				ady[i][j] = ady72[i][j];
			}
		}
	}
	
	//A�adir las aristas del chip
	if (N == 4)
		n = 4;
	else if (N == 8)
		n = 8;
	else if (N == 21)
		n = 24;
	else if (N == 40)
		n = 48;
	else if (N == 72)
		n = 121;
	
	// Iniciamos aristasDelChip
	aristasDelChip = new edge[n * 2];
	for (int i = 0; i < n; i++) {
		edge ed = edge();
		ed.pos1 = ady[i][0] - 1;
		ed.pos2 = ady[i][1] - 1;
		aristasDelChip[i] = ed;

		ed = edge();
		ed.pos1 = ady[i][1] - 1;
		ed.pos2 = ady[i][0] - 1;
		aristasDelChip[i + n] = ed;
	}
	
}

void IniciarQStates(){
	vector<int> a = puertas.a;
	vector<int> b = puertas.b;
	
	//Inicilizar L0 qstates

	qstatesTotal = new int[N]; // numero de puertas que ha de ejecutarse para cada qstate
	countTotal = new int[N];
	qstatesCount = new int[N];
	countCurrent = new int[N];
	

	for (int i = 0; i < N; i++) {
		qstatesCount[i] = 0;
		qstatesTotal[i] = 0;
		countCurrent[i] = 0;
		countTotal[i] = 0;
	}

	len = numTrabajos;
	numTareas = numTrabajos * ITERATIONS;

	//Datos para cromosoma
	
	puertasDelProblema = new operation[numTrabajos];

	// Iniciar puertasDelProblema qstatesTotal countTotal
	create_gnome(a, b, len, puertasDelProblema, ITERATIONS, qstatesTotal, countTotal);

	//Cota inferior de makespan
	int maxPS = 0;
	for (int i = 0; i < N; i++) {
		int auxc = 0;
		for (int j = 0; j < numTrabajos; j++) {
			if (a[j] == i || b[j] == i)
				auxc++;
		}
		if (auxc > maxPS)
			maxPS = auxc;
	}
	
	cotaInferior = maxPS * 3 * ITERATIONS + ITERATIONS;

	//Mejores saltos
	allAdy = new vector<int>[N];
	edgeAdy = new vector<int>*[N];
	bestHop = new vector<int>*[N];
	for (int i = 0; i < N; i++) {
		allAdy[i] = vector<int>();
		edgeAdy[i] = new vector<int>[N];
		bestHop[i] = new vector<int>[N];
		for (int j = 0; j < N; j++) {
			bestHop[i][j] = vector<int>();
			edgeAdy[i][j] = vector<int>();
		}
	}

	for (int i = 0; i < n; i++) {
		allAdy[ady[i][0] - 1].push_back(ady[i][1] - 1);
		allAdy[ady[i][1] - 1].push_back(ady[i][0] - 1);
	}


	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < allAdy[i].size(); k++) {
				if (distHop[i][j] - 1 == distHop[allAdy[i][k]][j])
					bestHop[i][j].push_back(allAdy[i][k]);
			}
		}
	}

	//Edges adyacentes
	for (int i = 0; i < n; i++) {
		edge ed = aristasDelChip[i];
		for (int j = 0; j < n; j++) {
			if (i != j) {
				edge ed2 = aristasDelChip[j];
				if (ed.pos1 == ed2.pos1 || ed.pos1 == ed2.pos2 || ed.pos2 == ed2.pos1 || ed.pos2 == ed2.pos2) {
					edgeAdy[ed.pos1][ed.pos2].push_back(j);
					edgeAdy[ed.pos1][ed.pos2].push_back(j + n);
					edgeAdy[ed.pos2][ed.pos1].push_back(j);
					edgeAdy[ed.pos2][ed.pos1].push_back(j + n);
				}
			}
		}
	}

	//PS adyacentes
	sharePS = new vector<int>[numTrabajos];
	for (int i = 0; i < numTrabajos; i++) {
		sharePS[i] = vector<int>();
	}

	for (int i = 0; i < numTrabajos; i++) {
		operation op = puertasDelProblema[i];
		for (int j = 0; j < numTrabajos; j++) {
			if (i != j) {
				operation op2 = puertasDelProblema[j];
				if (op.q1 == op2.q1 || op.q1 == op2.q2 || op.q2 == op2.q1 || op.q2 == op2.q2) {
					sharePS[i].push_back(j);
				}
			}
		}
	}
}

// create chromosome or string of genes 
void create_gnome(vector<int> q1, vector<int> q2, int len, operation* genOp, int iterations, int* qTotal, int* countCurrent)
{
	for (int i = 0; i < len; i++) {
		operation op = operation();
		op.type = "p-s";
		op.q1 = q1[i] - 1;
		op.q2 = q2[i] - 1;
		genOp[i] = op;
		qTotal[q1[i] - 1] = qTotal[q1[i] - 1]++; qTotal[q2[i] - 1] = qTotal[q2[i] - 1]++;
		countCurrent[q1[i] - 1] = countCurrent[q1[i] - 1]++; countCurrent[q2[i] - 1] = countCurrent[q2[i] - 1]++;
	}
}

//ReCalcular cabezas y tiempos de inicio
void BorraCabezasYColas(vector<operation>* operations) {

	for (int i = 0; i < operations->size(); i++) {
		(*operations)[i].head = -1;
		(*operations)[i].init = 0;
		(*operations)[i].tail = -1;
	}
}

// Inicia las estructuras
void InicializaEstructuras() {

	if ((METODO == 1) && (ELIMINARIGUALES > 0)) {
		makespanesPoblacion = new int[MAXIMOMAKESPAN];
		iteracionPoblacion = new int[MAXIMOMAKESPAN];
		for (int i = 0; i < MAXIMOMAKESPAN; i++) {
			makespanesPoblacion[i] = -1;
			iteracionPoblacion[i] = -1;
		}
	}


	int numeroinds;
	if (METODO == 0) numeroinds = TAMANOSOLUCIONESINICIALES;
	if (METODO == 2) numeroinds = TAMANOSOLUCIONESINICIALESTS;
	if ((METODO == 0) || (METODO == 2)) {
		SolucionesIniciales = new Individuo[numeroinds];
		for (int i = 0; i < numeroinds; i++) {
			IniciarIndividuo(SolucionesIniciales[i]);
		}
	}

	int numeroinds2;
	if (METODO == 0) numeroinds2 = TAMANOPOBLACION;
	if (METODO == 1) numeroinds2 = TAMANOPOBLACIONGENETICO;

	if ((METODO == 0) || (METODO == 1)) {
		Poblacion = new Individuo[numeroinds2];
		for (int i = 0; i < numeroinds2; i++) {
			IniciarIndividuo(Poblacion[i]);
		}
	}


	if (METODO == 1) {
		PoblacionAux = new Individuo[TAMANOPOBLACIONGENETICO];
		parejasAleatorias = new int[TAMANOPOBLACIONGENETICO];
		for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {
			IniciarIndividuo(PoblacionAux[i]);
		}
	}
	
	IniciarIndividuo(indActualBusquedaTabu);
	IniciarIndividuo(mejorIndividuoBusquedaTabu);
	IniciarIndividuo(individuoaux);
	IniciarIndividuo(individuoaux2);

	//Para depuracion y logs
	IniciarIndividuo(mejorIndividuoEncontrado);
	IniciarIndividuo(MejorIndividuoTotal);

	//Se inicializan las estructuras que utiliza la busqueda tabu
	vecinos = new int*[MAXVECINOS];
	arcosInvertidos = new int*[MAXVECINOS];
	for (int i = 0; i < MAXVECINOS; i++) {
		arcosInvertidos[i] = new int[5];
		vecinos[i] = new int[numTareas + 20];
	}
	estimaciones = new int[MAXVECINOS];

	//Utilizadas por el cruce JOX
	aleatorios = new int[numTrabajos];
	trabajosElegidos = new int[numTrabajos];
	if (INICIALES == 1) posElegidos = new int[N];

	MatrizDistanciasSolucionesIniciales = new int*[TAMANOSOLUCIONESINICIALES];
	for (int i = 0; i < TAMANOSOLUCIONESINICIALES; i++) MatrizDistanciasSolucionesIniciales[i] = new int[TAMANOSOLUCIONESINICIALES];

	MatrizDistancias = new int*[TAMANOPOBLACION];
	for (int i = 0; i < TAMANOPOBLACION; i++) MatrizDistancias[i] = new int[TAMANOPOBLACION];

	seleccionaux1 = 0;
	seleccionaux2 = 0;


}

// Iniciar individuo
void IniciarIndividuo(Individuo &ind) {
	ind.vectorSecuenciaOperaciones = new int[numTareas];
	ind.vectorPosicionOperaciones = new int[numTareas];
//	if (INICIALES == 1) ind.iniciales = new int[N];
	ind.vectorSecuenciaEdges = new int[numTareas];
	ind.qbitsBase = new int[N];
	ind.qstatesBase = new int[N];
}

int existeElem (int* v, int longi, int num)
{
	num = num ;
    int i, enc=0;
    for(i=0; i<longi && !enc; i++)
    {
        if(v[i]==num)
            enc=1;
    }
    return enc;
}

// COSAS TFM
void generarVector (int* v)
{
    int i, j, num, longi=N;
    srand(time(NULL));
    for (i=0; i<longi; i++)
    {
        while(existeElem(v, longi, num = rand() % 4 + 1 ));
        v[i]= num;
    }
}

void restarUno (int* v, int longi)
{
    int i;
    for(i=0; i<longi; i++)
    	v[i] = v[i] - 1;
}

void escribirVec (int* v, int longi)
{
    int i;
    longi=N;
    for(i=0; i<longi; i++)
        printf("%d ", v[i]);
}


//IMPORTANTE
void CrearPoblacionInicialGenetico() {

	// Creamos individuos aleatorios
	for (int i = 0; i < TAMANOPOBLACIONGENETICO; i++) {

		CrearVectoresIniciales(Poblacion[i]);


		//Para cada solucion inicial, la planificamos, calculamos las colas, generamos caminos criticos,
		// hacemos busqueda tabu y construimos los grafos
		PlanificacionSemiactiva(Poblacion[i], itActual, itFinal);

	}
}

//Crea poblacion inicial
void CrearVectoresIniciales(Individuo &ind) {

	//Crear vectores de secuenciamiento aleatorios
	//Vectores de ps gates
	for (int j = itActual * numTrabajos; j < numTareas; j++)
		ind.vectorSecuenciaOperaciones[j] = j;
	
	for (int i = itActual; i < ITERATIONS; i++) {
		for (int j = i * numTrabajos; j < i * numTrabajos + numTrabajos; j++) {
			int nuevaPosicionAleatoria = Aleatorio(i * numTrabajos, i * numTrabajos + numTrabajos - 1);
			int posicionbackup = ind.vectorSecuenciaOperaciones[j];
			ind.vectorSecuenciaOperaciones[j] = ind.vectorSecuenciaOperaciones[nuevaPosicionAleatoria];
			ind.vectorSecuenciaOperaciones[nuevaPosicionAleatoria] = posicionbackup;
		}
	}

	//Vectores posicion de la secuencias de puertas ps
	for (int j = itActual * numTrabajos; j < numTareas; j++) {
		ind.vectorPosicionOperaciones[ind.vectorSecuenciaOperaciones[j]] = j;
	}

	//iniciar una solucion incial regular siguiendo un heuristico
	if (primero >= 0) 
	{

		//Inicilizar L0 qstates
		primero--;

		int* qbits = new int[N];
		int* qstates = new int[N];
		if (itActual == 0) 
		{
			for (int j = 0; j < N; j++) 
			{
//				qbits[j] = j;
//				qstates[j] = j;
				qbits[j] = vectorInicios[j];
				qstates[vectorInicios[j]] = j;
			}
		}
		else 
		{
			for (int j = 0; j < N; j++) 
			{
				qbits[j] = ind.qbitsBase[j];
				qstates[j] = ind.qstatesBase[j];
			}
		}
		for (int j = itActual * numTrabajos; j < itFinal*numTrabajos; j++) 
		{
			operation op = puertasDelProblema[ind.vectorSecuenciaOperaciones[j] % numTrabajos];//get op
			//encontra un edge adyacente a op.q1 a partir de ady y genedge
			int i = 0;
			int ps[] = { qstates[op.q1],qstates[op.q2] };
			int b = bestHop[ps[i % 2]][ps[(i + 1) % 2]][0];
			while (b != ps[(i + 1) % 2])
			{
				int a = ps[i % 2];
				//swap qstates
				int q = qbits[a];
				qbits[a] = qbits[b];
				qbits[b] = q;

				qstates[qbits[a]] = a;
				qstates[qbits[b]] = b;

				ps[i % 2] = b;
				i++;
				b = bestHop[ps[i % 2]][ps[(i + 1) % 2]][0];
			}
			i = 0;
			while (i < n && !(aristasDelChip[i].pos1 == ps[0]
				&& aristasDelChip[i].pos2 == ps[1] || aristasDelChip[i].pos1 == ps[1]
				&& aristasDelChip[i].pos2 == ps[0]))
				i++;
			ind.vectorSecuenciaEdges[j] = i;
		}
	}
	else {
		for (int j = itActual * numTrabajos; j < numTareas; j++) 
		{
			int nuevaPosicionAleatoria = Aleatorio(0, n * 2 - 1);
			ind.vectorSecuenciaEdges[j] = nuevaPosicionAleatoria;
		}
	}
}

//DEvuelve un numero aleatorio
int Aleatorio(int minimo, int maximo) {
	return rand() % (maximo - minimo + 1) + minimo;
}

//IMPORTANTE es crear la solucion
void PlanificacionSemiactiva(Individuo &ind, int itIni, int itFin) {
	//Inicializar Schedule
	vector<operation> operations = vector<operation>();
	operations.reserve(numTareas * 10);
	//Planificacion
	PlanificacionSemiactiva(ind, &operations, itIni, itFin);

	if (guardar) {
		ind.operations = mejoresInicios.size();
		mejoresInicios.push_back(operations);
	}
	else if (guardar) {
		ind.operations = guardados;
		guardados++;
		mejoresInicios[ind.operations] = operations;
	}

	//Free
	vector<operation>().swap(operations);
}

void PlanificacionSemiactiva(Individuo &ind, vector<operation> *operations, int itIni, int itFin) {

	//Inicializar Schedule
	if (itIni > 0 && SOLOBL) {
		*operations = mejoresInicios[ind.operations];
	}
	//else ind.numerodeswaps = 0;
	GeneraSchedule(operations, ind, itIni, itFin);

	//Makespan
	int makespan_f = 0;
	//TotalHead(operations);
	for (int i = 0; i < N; i++) {
		int m;
		m = (*operations)[(*operations)[i].pre1].init + (*operations)[(*operations)[i].pre1].time;
		if (m > makespan_f)
			makespan_f = m;
	}

	ind.makespan = makespan_f + 500 * (ITERATIONS - itFin);

}

//genera una planificacion
void GeneraSchedule(std::vector<operation> *operations, Individuo & ind, int itIni, int itFin)
{
	//Inicilizar L0 qstates
	int* qbits = new int[N];
	int* qstates = new int[N];
	int* qCount = new int[N];
//	cout << endl;
	for (int i = 0; i < N; i++) {
		operation op = operation();
		if (itIni == 0) {
			op.q1 = i;
			op.q2 = i;
			op.qb1 = qstates[i];
			op.qb2 = qstates[i];
			op.orden = i;
			op.tail = -1;
			op.head = 0;
			op.init = 0;
			op.time = 0;
			op.suc1 = i;
			op.pre1 = i;
			op.suc2 = i;
			op.pre2 = i;
			operations->push_back(op);
			qCount[i] = 0;
			countCurrent[i] = 0;
			if (INICIALES == 1) {
				qbits[ind.iniciales[i]] = i;
				qstates[i] = ind.iniciales[i];
			}
			// Inicializacion de qbits Trabajo a implementar por mi (TFM)
			else {
//				qbits[i] = i;
//				qstates[i] = i;
				qstates[vectorInicios[i]] = i;
				qbits[i] = vectorInicios[i];
//				cout << "El QBit: " << qbits[i] << " Tiene QState: " << qstates[i] << endl;
			}
		}
		else {
			qbits[i] = ind.qbitsBase[i];
			qstates[i] = ind.qstatesBase[i];
			qCount[i] = qstatesTotal[i] * itIni;
			countCurrent[i] = qstatesTotal[i] * itIni;
		}
	}

	//reordenar si es necesario
	//regen(&ind, countTotal, countCurrent, regenOp, len*ITERATIONS);

	//while fin
	int i = itIni * numTrabajos;
	operation op; edge ed;
	while (i < len*itFin) {
		// Modificacion a posteriori de la inicializacion (TFM)
		//Get first op and edge
		op = puertasDelProblema[ind.vectorSecuenciaOperaciones[i] % numTrabajos];
		op.posCrh = i;
		ed = aristasDelChip[ind.vectorSecuenciaEdges[i]];
		//Insert Swaps
		insertSWAP(ed, &op, operations, qbits, N, qstates);
		//InsertOp
		insertPS(&op, operations, qbits, qstates, N, qstatesTotal, qCount);

		if (i == len - 1 && guardar) {
			//Preparar mejor
			for (int j = 0; j < N; j++) {
				ind.qbitsBase[j] = qbits[j];
				ind.qstatesBase[j] = qstates[j];
			}
		}

		i++;
	}

	if (SELECCION == 4) {
		ind.qstatesBaseSuma = 0;
		for (int j = 0; j < N; j++)
			ind.qstatesBaseSuma += qstates[j] * pow(N, j);
	}


	//free
	delete qbits;
	delete qstates;
	delete qCount;
}

//Añade swaps si son necesarios para poder colocar una puerta PS
void insertSWAP(edge ed, operation *op, vector<operation>* operations, int *qbits, int len, int *qstates) {
	//puntos de inicio
	int t1 = 0;
	int t2 = 0;
	//qstate ubicado en qbit
	int i1 = 0;
	int i2 = 0;
	//Necesita swap
	bool needSwap = true;

	//insertar hasta ue no sea necesario
	while (needSwap) {
		//Localizar qstates
		LocalizarQstates(qbits, qstates, operations, op, &t1, &i1, &t2, &i2);
		//Comprobar si es posible aÃ±adir operacion
		if (i1 == ed.pos1 && i2 == ed.pos2 || i1 == ed.pos2 && i2 == ed.pos1) {
			//No necesita swap
			needSwap = false;
		}
		else {
			//Coloca swaps de forma que el camino mas largo sea el minimo posible
			if (max(distHop[i2][ed.pos1], distHop[i1][ed.pos2]) <
				max(distHop[i1][ed.pos1], distHop[i2][ed.pos2])) {
				//	S â† InsertSwap(menor distancia)
				int i = i2, p = ed.pos1, o = i1, po = ed.pos2;
				if (i1 == ed.pos2) {
					i = i2, p = ed.pos1, o = i1, po = ed.pos2;
				}

				else if (distHop[i2][ed.pos1] < distHop[i1][ed.pos2])
					i = i1, p = ed.pos2, o = i2, po = ed.pos1;
				//Siguiente salto
				//cout << i << " " << p << "a\n";
				int si = bestHop[i][p][operations->size() % bestHop[i][p].size()];
				if (si == o && o != po) {
					//cout << o << " " << po << "e\n";
					i = o; si = bestHop[i][po][operations->size() % bestHop[i][po].size()];
				}
				else if (si == o && o != p) {
					//cout << o << " " << p << "f\n";
					i = o; si = bestHop[i][p][operations->size() % bestHop[i][p].size()];
				}
				//add swap
				addSwap(i, si, ed, operations, qbits, qstates);
			}
			else {
				//	S â† InsertSwap(menor distancia)
				int i = i1, o = i2, p = ed.pos1, po = ed.pos2;
				if (i2 == ed.pos2) {
					i = i1, o = i2, p = ed.pos1, po = ed.pos2;
				}
				else if (distHop[i1][ed.pos1] < distHop[i2][ed.pos2]) {
					i = i2, p = ed.pos2, o = i1, po = ed.pos1;
				}
				//Siguiente salto
				//cout << i << " " << p << "d\n";
				int si = bestHop[i][p][operations->size() % bestHop[i][p].size()];
				if (si == o && o != po) {
					//cout << o << " " << po << "e\n";
					i = o; si = bestHop[i][po][operations->size() % bestHop[i][po].size()];
				}
				else if (si == o && o != p) {
					//cout << o << " " << p << "f\n";
					i = o; si = bestHop[i][p][operations->size() % bestHop[i][p].size()];
				}
				//add swap
				addSwap(i, si, ed, operations, qbits, qstates);
			}
		}
	}
}

//Anade una puerta PS
void insertPS(operation *op, vector<operation>* operations, int *qbits, int *qstates, int len, int* qstatesTotal, int* qstatesCount) {
	//puntos de inicio
	int t1 = 0;
	int t2 = 0;
	//qstate ubicado en qbit
	int i1 = 0;
	int i2 = 0;
	//Localizar qstates
	LocalizarQstates(qbits, qstates, operations, op, &t1, &i1, &t2, &i2);

	//Init op
	int t = t2;
	if (t1 > t2)
		t = t1;


	//Â¿hay hueco?
	(*op).init = t;
	(*op).time = timePS[i1][i2];
	(*op).tail = -1;
	int tf = (*op).time + (*op).init;
	(*op).qb1 = qstates[(*op).q1];
	(*op).qb2 = qstates[(*op).q2];
	//Sucesores y predecesores
	(*op).suc1 = (*op).q1;
	(*op).suc2 = (*op).q2;
	(*op).pre1 = (*operations)[(*op).q1].pre1;
	(*op).pre2 = (*operations)[(*op).q2].pre1;

	(*op).head = -1;
	//Add op
	int pos = (int)(*operations).size();
	(*op).orden = pos;
	(*op).pos = pos;
	operations->push_back(*op);
	qstatesCount[(*op).q1]++;
	qstatesCount[(*op).q2]++;

	//sucesores
	if ((*operations)[(*op).pre1].suc1 < N || (*operations)[(*op).pre1].type == "mix") {
		(*operations)[(*op).pre1].suc1 = pos;
		if ((*operations)[(*op).pre1].type == "mix" || (*op).pre1 < N)
			(*operations)[(*op).pre1].suc2 = pos;
	}
	else
		(*operations)[(*op).pre1].suc2 = pos;

	if ((*operations)[(*op).pre2].suc1 < N || (*operations)[(*op).pre2].type == "mix") {
		(*operations)[(*op).pre2].suc1 = pos;
		if ((*operations)[(*op).pre2].type == "mix" || (*op).pre2 < N)
			(*operations)[(*op).pre2].suc2 = pos;
	}
	else
		(*operations)[(*op).pre2].suc2 = pos;

	//predecesor de final
	(*operations)[(*op).q1].pre1 = pos;
	(*operations)[(*op).q2].pre1 = pos;
	(*operations)[(*op).q1].pre2 = pos;
	(*operations)[(*op).q2].pre2 = pos;

	//Mix?
	addMix(qstatesCount, (*op).q1, qstatesTotal, tf, operations, i1, qstates);
	addMix(qstatesCount, (*op).q2, qstatesTotal, tf, operations, i2, qstates);
}


//Copia un individuo
void CopiarIndividuo(Individuo &ind1, Individuo &ind2, int nivel) {

	if (nivel == 0) {
		for (int j = 0; j < numTareas; j++) ind2.vectorSecuenciaOperaciones[j] = ind1.vectorSecuenciaOperaciones[j];
	}

	if (nivel >= 1) {
		ind2.makespan = ind1.makespan;
		for (int j = 0; j < numTrabajos*ITERATIONS; j++) 
		{
			ind2.vectorSecuenciaOperaciones[j] = ind1.vectorSecuenciaOperaciones[j];
			ind2.vectorPosicionOperaciones[j] = ind1.vectorPosicionOperaciones[j];
			ind2.vectorSecuenciaEdges[j] = ind1.vectorSecuenciaEdges[j];
		}

		for (int j = 0; j < N; j++) 
		{
			ind2.qbitsBase[j] = ind1.qbitsBase[j];
			ind2.qstatesBase[j] = ind1.qstatesBase[j];
		}

		ind2.qstatesBaseSuma = ind1.qstatesBaseSuma;
		ind2.operations = ind1.operations;
		if (INICIALES == 1) {
			for (int j = 0; j < N; j++)
				ind2.iniciales[j] = ind1.iniciales[j];
		}
	}
}






//Localiza la posicion de los qstates
void LocalizarQstates(int * qbits, int * qstates, vector<operation> *operations, operation *op, int *t1, int *i1, int *t2, int *i2)
{

	*i1 = qstates[(*op).q1];
	operation lastOP = (*operations)[(*operations)[(*op).q1].pre1];
	*t1 = lastOP.init + lastOP.time;

	*i2 = qstates[(*op).q2];
	lastOP = (*operations)[(*operations)[(*op).q2].pre1];
	*t2 = lastOP.init + lastOP.time;
}

//Añade un swap entre dos posiciones y cambia la posicion de los qstates
void addSwap(int i, int si, edge &ed, std::vector<operation> *operations, int *qbits, int *qstates)
{
	//Get qstates
	int qi = qbits[i];
	int qsi = qbits[si];

	//add swap
	operation swap;
	swap.type = "swap";
	swap.time = TIME_SWAP;
	swap.q1 = qi;
	swap.q2 = qsi;
	swap.pos = -1;
	swap.qb1 = qstates[swap.q1];
	swap.qb2 = qstates[swap.q2];

	//t inicio
	int t, ts = 0, ti = 0;
	ti = (*operations)[(*operations)[qi].pre1].init + (*operations)[(*operations)[qi].pre1].time;
	ts = (*operations)[(*operations)[qsi].pre1].init + (*operations)[(*operations)[qsi].pre1].time;
	if (ti > ts)
		t = ti;
	else
		t = ts;


	swap.init = t;

	//Sucesores y predecesores
	swap.suc1 = swap.q1;
	swap.suc2 = swap.q2;
	swap.pre1 = (*operations)[swap.q1].pre1;
	swap.pre2 = (*operations)[swap.q2].pre1;

	swap.tail = -1;
	swap.head = -1;

	//revert mix
	bool mix1 = false, mix2 = false;
	int mixi1 = 0, mixi2 = 0;

	if ((*operations)[(*operations)[qi].pre1].type == "mix" && (*operations)[(*operations)[qi].pre1].init + (*operations)[(*operations)[qi].pre1].time == t) {
		if ((*operations)[(*operations)[qsi].pre1].init + (*operations)[(*operations)[qsi].pre1].time <= t - 1
			|| ((*operations)[(*operations)[qsi].pre1].init + (*operations)[(*operations)[qsi].pre1].time == t && (*operations)[(*operations)[qsi].pre1].type == "mix")) {
			mix1 = true;
			//(*operations)[(*operations)[qi].pre1].orden = -1;
			mixi1 = (*operations)[qi].pre1;
			swap.pre1 = (*operations)[(*operations)[qi].pre1].pre1;
			swap.init = t - 1;
		}
	}


	else if ((*operations)[(*operations)[qsi].pre1].type == "mix" && (*operations)[(*operations)[qsi].pre1].init + (*operations)[(*operations)[qsi].pre1].time == t) {
		if ((*operations)[(*operations)[qi].pre1].init + (*operations)[(*operations)[qi].pre1].time <= t - 1
			|| ((*operations)[(*operations)[qi].pre1].init + (*operations)[(*operations)[qi].pre1].time == t && (*operations)[(*operations)[qi].pre1].type == "mix")) {
			mix2 = true;
			//(*operations)[(*operations)[qsi].pre1].orden = -1;
			mixi2 = (*operations)[qsi].pre1;
			swap.pre2 = (*operations)[(*operations)[qsi].pre1].pre1;
			swap.init = t - 1;
		}
	}


	//add swap
	int pos = (int)(*operations).size();
	swap.orden = pos;
	swap.pos = pos;
	operations->push_back(swap);

	//sucesores
	if ((*operations)[swap.pre1].suc1 < N || (*operations)[swap.pre1].type == "mix" || (mix1 && (*operations)[swap.pre1].suc1 == mixi1)) {
		(*operations)[swap.pre1].suc1 = pos;
		if (swap.pre1 < N || (*operations)[swap.pre1].type == "mix")
			(*operations)[swap.pre1].suc2 = pos;
	}
	else
		(*operations)[swap.pre1].suc2 = pos;

	if ((*operations)[swap.pre2].suc1 < N || (*operations)[swap.pre2].type == "mix" || (mix2 && (*operations)[swap.pre2].suc1 == mixi2)) {
		(*operations)[swap.pre2].suc1 = pos;
		if (swap.pre2 < N || (*operations)[swap.pre2].type == "mix")
			(*operations)[swap.pre2].suc2 = pos;
	}
	else
		(*operations)[swap.pre2].suc2 = pos;


	//predecesor de final
	(*operations)[swap.q1].pre1 = pos;
	(*operations)[swap.q2].pre1 = pos;
	(*operations)[swap.q1].pre2 = pos;
	(*operations)[swap.q2].pre2 = pos;

	//swap qstates
	int q = qbits[i];
	qbits[i] = qbits[si];
	qbits[si] = q;

	qstates[qbits[i]] = i;
	qstates[qbits[si]] = si;

	//add mix
	if (mix1) {
		addMix(qstatesCount, qi, qstatesTotal, swap.init + TIME_SWAP, operations, 0, qstates);
	}
	if (mix2) {
		addMix(qstatesCount, qsi, qstatesTotal, swap.init + TIME_SWAP, operations, 0, qstates);
	}
}

//Añade un mix si es necesario
void addMix(int * qstatesCount, int qmix, int * qstatesTotal, int ti, std::vector<operation> *operations, int i, int* qstates)
{
	if (qstatesTotal[qmix] != 0 && qstatesCount[qmix] % qstatesTotal[qmix] == 0) {
		operation mix = operation();
		mix.q1 = qmix;
		mix.qb1 = qstates[mix.q1];
		mix.qb2 = qstates[mix.q1];
		mix.type = "mix";
		mix.time = TIME_MIX;
		mix.tail = -1;
		mix.pre1 = (*operations)[mix.q1].pre1;
		mix.pre2 = (*operations)[mix.q1].pre1;
		mix.suc1 = mix.q1;
		mix.suc2 = mix.q1;
		mix.head = -1;
		int t = ti;


		mix.init = t;
		int pos = (int)(*operations).size();
		mix.orden = pos;
		mix.pos = pos;
		operations->push_back(mix);
		//Sucesores de predecesor
		if ((*operations)[mix.pre1].suc1 < N)
			(*operations)[mix.pre1].suc1 = pos;
		else
			(*operations)[mix.pre1].suc2 = pos;

		//predecesor de final
		(*operations)[qmix].pre1 = pos;
		(*operations)[qmix].pre2 = pos;
	}
}

//Calcular cola
int Cola(int op, vector<operation>* operations) {
	int opAct = op;
	vector<int> opPendant = vector<int>();
	opPendant.push_back(opAct);

	while (!opPendant.empty()) {
		opAct = opPendant.back();

		if ((*operations)[opAct].suc1 < N && (*operations)[opAct].suc2 < N) {
			(*operations)[opAct].tail = 0;
			opPendant.pop_back();
		}
		else if ((*operations)[opAct].tail <= -1) {
			if ((*operations)[opAct].suc1 >= N && (*operations)[opAct].suc2 >= N && (*operations)[opAct].type != "mix" && opAct >= N) {

				if ((*operations)[(*operations)[opAct].suc1].tail == -1) {
					(*operations)[(*operations)[opAct].suc1].tail = -2;
					opPendant.push_back((*operations)[opAct].suc1);
				}
				else if ((*operations)[(*operations)[opAct].suc2].tail == -1) {
					(*operations)[(*operations)[opAct].suc2].tail = -2;
					opPendant.push_back((*operations)[opAct].suc2);
				}
				else {
					int cola = (*operations)[(*operations)[opAct].suc1].tail;
					int time = (*operations)[(*operations)[opAct].suc1].time;
					int cola2 = (*operations)[(*operations)[opAct].suc2].tail;
					int time2 = (*operations)[(*operations)[opAct].suc2].time;
					(*operations)[opAct].tail = max(cola + time, cola2 + time2);
					opPendant.pop_back();
				}
			}
			else if ((*operations)[opAct].suc1 >= N) {
				if ((*operations)[(*operations)[opAct].suc1].tail == -1) {
					(*operations)[(*operations)[opAct].suc1].tail = -2;
					opPendant.push_back((*operations)[opAct].suc1);
				}
				else {
					int cola = (*operations)[(*operations)[opAct].suc1].tail;
					int time = (*operations)[(*operations)[opAct].suc1].time;
					(*operations)[opAct].tail = cola + time;
					opPendant.pop_back();
				}
			}
			else if ((*operations)[opAct].suc2 >= N && (*operations)[opAct].type != "mix" && opAct >= N) {
				if ((*operations)[(*operations)[opAct].suc2].tail == -1) {
					(*operations)[(*operations)[opAct].suc2].tail = -2;
					opPendant.push_back((*operations)[opAct].suc2);
				}
				else {
					int cola = (*operations)[(*operations)[opAct].suc2].tail;
					int time = (*operations)[(*operations)[opAct].suc2].time;
					(*operations)[opAct].tail = cola + time;
					opPendant.pop_back();
				}
			}
		}
		else {
			opPendant.pop_back();
		}
	}
	vector<int>().swap(opPendant);

	return (*operations)[op].tail;
}


//Calcular cabeza
int Cabeza(int op, vector<operation>* operations) {
	vector<int> opPendant = vector<int>();
	int opAct = op;
	opPendant.push_back(opAct);

	while (!opPendant.empty()) {
		opAct = opPendant.back();


		if ((*operations)[opAct].pre1 < N && (*operations)[opAct].pre2 < N) {
			(*operations)[opAct].head = 0;
			(*operations)[opAct].init = 0;
			opPendant.pop_back();
		}
		else if ((*operations)[opAct].head <= -1) {
			if ((*operations)[opAct].pre2 >= N && (*operations)[opAct].pre1 >= N && (*operations)[opAct].type != "mix") {

				if ((*operations)[(*operations)[opAct].pre1].head == -1) {
					(*operations)[(*operations)[opAct].pre1].head = -2;
					opPendant.push_back((*operations)[opAct].pre1);
				}
				else if ((*operations)[(*operations)[opAct].pre2].head == -1) {
					(*operations)[(*operations)[opAct].pre2].head = -2;
					opPendant.push_back((*operations)[opAct].pre2);
				}
				else {
					int head = (*operations)[(*operations)[opAct].pre1].head;
					int time = (*operations)[(*operations)[opAct].pre1].time;
					int init = (*operations)[(*operations)[opAct].pre1].init;
					int head2 = (*operations)[(*operations)[opAct].pre2].head;
					int time2 = (*operations)[(*operations)[opAct].pre2].time;
					int init2 = (*operations)[(*operations)[opAct].pre2].init;
					(*operations)[opAct].head = max(head + time, head2 + time2);
					(*operations)[opAct].init = max(init + time, init2 + time2);
					opPendant.pop_back();
				}
			}
			else if ((*operations)[opAct].pre2 >= N && (*operations)[opAct].type != "mix" && opAct >= N) {
				if ((*operations)[(*operations)[opAct].pre2].head == -1) {
					(*operations)[(*operations)[opAct].pre2].head = -2;
					opPendant.push_back((*operations)[opAct].pre2);
				}
				else {
					int head = (*operations)[(*operations)[opAct].pre2].head;
					int time = (*operations)[(*operations)[opAct].pre2].time;
					int init = (*operations)[(*operations)[opAct].pre2].init;
					(*operations)[opAct].head = head + time;
					(*operations)[opAct].init = init + time;
					opPendant.pop_back();
				}
			}
			else if ((*operations)[opAct].pre1 >= N) {
				if ((*operations)[(*operations)[opAct].pre1].head == -1) {
					(*operations)[(*operations)[opAct].pre1].head = -2;
					opPendant.push_back((*operations)[opAct].pre1);
				}
				else {

					int head = (*operations)[(*operations)[opAct].pre1].head;
					int time = (*operations)[(*operations)[opAct].pre1].time;
					int init = (*operations)[(*operations)[opAct].pre1].init;
					(*operations)[opAct].head = head + time;
					(*operations)[opAct].init = init + time;

					opPendant.pop_back();
				}
			}
		}
		else {
			opPendant.pop_back();
		}
	}

	vector<int>().swap(opPendant);

	return (*operations)[op].head;
}



//Calcular el siguiente salto en el camino critico
int SiguienteCaminoCritico(int ope, vector<operation>* operations, int makespan) {
	int c = -1; //Head(ope, operations);
	Cola(ope, operations);
	operation op = (*operations)[ope];
	if (op.pre1 < N && op.pre2 < N)
		return c;
	else if (op.orden < N)
		return op.pre1;
	else {
		if (op.pre2 >= N && op.pre1 >= N) {
			int h = op.tail;
			int m1 = h + op.time + (*operations)[op.pre1].head + (*operations)[op.pre1].time;
			int m2 = h + op.time + (*operations)[op.pre2].head + (*operations)[op.pre2].time;
			if (m1 == m2 == makespan) {
				if ((*operations)[op.pre2].type == "p-s")
					//if (ope % 2 == 0)
					c = op.pre2;
				else
					c = op.pre1;
			}
			else if (m1 == makespan)
				c = op.pre1;
			else if (m2 == makespan)
				c = op.pre2;
			else
				exit(7);
		}
		else if (op.pre2 >= N && op.tail + op.time + (*operations)[op.pre2].head + (*operations)[op.pre2].time == makespan) {
			c = op.pre2;
		}
		else if (op.tail + op.time + (*operations)[op.pre1].head + (*operations)[op.pre1].time == makespan) {
			c = op.pre1;
		}
		else
			exit(8);
	}
	return c;
}

//Escribe la planificacion
void PrintSchedule2(std::vector<operation> &operations)
{
	cout << "\n";
	cout << "id; type; qag1; qag2; qloc1; qloc2; st; et; just1; just2" << "\n";
	//Makespan
	int j = 0;
	for (int i = 0; i < operations.size(); i++) {
		if (i < N)
			cout << j << ";" << operations[i].type << ";" <<
			operations[i].q1 << ";" << 0 << ";" <<
			operations[i].qb1 << ";" << 0 << ";" <<
			0 << ";" << 0 << ";" <<
			i << ";" << -1 << "\n";
		else if (operations[i].type == "MIX")
			cout << j << ";" << operations[i].type << ";" <<
			operations[i].q1 << ";" << 0 << ";" <<
			operations[i].qb1 << ";" << 0 << ";" <<
			operations[i].init << ";" << operations[i].time << ";" <<
			operations[i].pre1 << ";" << -1 << "\n";
		else if (operations[i].type != "")
			cout << j << ";" << operations[i].type << ";" <<
			operations[i].q1 << ";" << operations[i].q2 << ";" <<
			operations[i].qb1 << ";" << operations[i].qb2 << ";" <<
			operations[i].init << ";" << operations[i].time << ";" <<
			operations[i].pre1 << ";" << operations[i].pre2 << "\n";
		if (operations[i].type != "") j++;
	}
}

//Escribe la planificacion por qstates
void PrintSchedule(std::vector<operation> &operations)
{
	cout << "\n";
	//cout << "Caminos:\n";
	//Makespan
	for (int i = 0; i < N; i++) {
		if (operations.at(i).pre1 == i)
			cout << "Makespan " << i + 1 << ": " << 0 << "\n";
		else {
			int cola = Cola(i, &operations);
			int head = Cabeza(i, &operations);
			int makespan = operations.at(operations.at(i).pre1).init + operations.at(operations.at(i).pre1).time;
			cout << "Makespan " << i + 1 << ": " << makespan << "\n";
			//cout << "Cola: " << cola << "\n";
			//cout << "Head: " << head << "\n";
			int j = i;
			do {
				operation opj = operations.at(j);
				if (operations.at(opj.suc1).q1 == i || (operations.at(opj.suc1).type != "MIX" && operations.at(opj.suc1).q2 == i))
					j = opj.suc1;
				else if (operations.at(opj.suc2).q1 == i || (operations.at(opj.suc2).type != "MIX" && operations.at(opj.suc2).q2 == i))
					j = opj.suc2;
				else { cout << "F"; exit(3); }
				opj = operations.at(j);
				if (j >= N) {
					if (opj.type != "MIX")
						cout << opj.type.c_str() << " q1: " << opj.q1 + 1 << " q2: " << opj.q2 + 1 << " qb1: " << opj.qb1 + 1 << " qb2: " << opj.qb2 + 1 << " init: " << opj.init << "  time: " << opj.time << "; \n";
					else
						cout << opj.type.c_str() << ": " << opj.q1 + 1 << " qb1: " << opj.qb1 + 1 << " init: " << opj.init << "  time: " << opj.time << "; \n";
				}
			} while (j >= N);
			cout << "\n";
		}
	}
}


//Que CruceJOX deje dos individuos, uno en PoblacionAux[ind1] y otro en PoblacionAux[ind2]
void Cruce(int ind1, int ind2) {
	//int c = (numeroiteraciones%2)/1;
	if (TIPOCRUCE == 0) {    //Cruce JOX
		//Iteracion
		int it = itActual;
		while (it < itFinal) {
			//Indices
			int indice1 = it * numTrabajos;
			int indice2 = it * numTrabajos;
			int limite = it * numTrabajos + numTrabajos;

			//Elegidos
			int numTrab = Aleatorio(0, numTrabajos - 1);
			for (int i = 0; i < numTrabajos; i++) aleatorios[i] = i;
			for (int ind = 0; ind < numTrabajos; ind++) {
				trabajosElegidos[ind] = 0;
				int aux = Aleatorio(0, numTrabajos - 1);
				int temp = aleatorios[ind];
				aleatorios[ind] = aleatorios[aux];
				aleatorios[aux] = temp;
			}
			for (int i = 1; i <= numTrab; i++) {
				trabajosElegidos[aleatorios[i - 1]] = 1;
			}

			//Cruce
			while ((indice1 < limite) && (indice2 < limite)) {
				int cambiar = 1;
				//cout << ind1 << " " << ind2 << " " << indice1 << " " << PoblacionAux[ind1].vectorSecuencia[indice1] << " ";
				if (trabajosElegidos[PoblacionAux[ind1].vectorSecuenciaOperaciones[indice1] - numTrabajos * it] == 1) {
					indice1++;
					cambiar = 0;
				}
				else if (trabajosElegidos[PoblacionAux[ind2].vectorSecuenciaOperaciones[indice2] - numTrabajos * it] == 1) {
					indice2++;
					cambiar = 0;
				}
				if (cambiar == 1) {
					int aux = PoblacionAux[ind1].vectorSecuenciaOperaciones[indice1];
					int aux2 = PoblacionAux[ind1].vectorSecuenciaEdges[indice1];
					PoblacionAux[ind1].vectorSecuenciaOperaciones[indice1] = PoblacionAux[ind2].vectorSecuenciaOperaciones[indice2];
					PoblacionAux[ind1].vectorSecuenciaEdges[indice1] = PoblacionAux[ind2].vectorSecuenciaEdges[indice2];
					PoblacionAux[ind1].vectorPosicionOperaciones[PoblacionAux[ind2].vectorSecuenciaOperaciones[indice2]] = indice1;
					PoblacionAux[ind2].vectorSecuenciaOperaciones[indice2] = aux;
					PoblacionAux[ind2].vectorSecuenciaEdges[indice2] = aux2;
					PoblacionAux[ind2].vectorPosicionOperaciones[aux] = indice2;
					indice1++;
					indice2++;
				}
			}


			it++;
		}

		if (INICIALES == 1) {
			int numTrab = Aleatorio(1, N - 1);
			for (int i = 0; i < N; i++) aleatorios[i] = i;
			for (int ind = 0; ind < N; ind++) {
				posElegidos[ind] = 0;
				int aux = Aleatorio(0, N - 1);
				int temp = aleatorios[ind];
				aleatorios[ind] = aleatorios[aux];
				aleatorios[aux] = temp;
			}
			for (int i = 1; i <= numTrab; i++) {
				posElegidos[aleatorios[i - 1]] = 1;
			}
			//==========
			//  for(int i=0; i<numTrabajos; i++) trabajosElegidos[i]=0;
			//  trabajosElegidos[Aleatorio(0, numTrabajos-1)]=1;
			//==========

			int itCruce = 0;

			int indice1 = 0;
			int indice2 = 0;
			int limite = N;
			while ((indice1 < limite) && (indice2 < limite)) {
				int cambiar = 1;
				//cout << ind1 << " " << ind2 << " " << indice1 << " " << PoblacionAux[ind1].vectorSecuencia[indice1] << " ";
				if (posElegidos[PoblacionAux[ind1].iniciales[indice1]] == 1) {
					indice1++;
					cambiar = 0;
				}
				else if (posElegidos[PoblacionAux[ind2].iniciales[indice2]] == 1) {
					indice2++;
					cambiar = 0;
				}
				if (cambiar == 1) {
					int aux = PoblacionAux[ind1].iniciales[indice1];
					PoblacionAux[ind1].iniciales[indice1] = PoblacionAux[ind2].iniciales[indice2];
					PoblacionAux[ind2].iniciales[indice2] = aux;
					//cout << PoblacionAux[ind1].iniciales[indice1] << " " << PoblacionAux[ind2].iniciales[indice2] << "\n";
					indice1++;
					indice2++;
				}
			}

		}

	}

	if (TIPOCRUCE == 1) {    //Cruce NWOX

		//Iteracion
		int it = itActual;
		//Indices
		int indice1 = it * numTrabajos;
		int limite;
		if (itFinal == ITERATIONS && itActual == 0)
			limite = numTareas;
		else
			limite = it * numTrabajos + numTrabajos;

		//Elegimos dos posiciones aleatorias posmin y posmax
		int posmin = Aleatorio(indice1, limite - 1);
		int posmax = Aleatorio(indice1, limite - 1);
		if (posmin > posmax) {
			int aux = posmax;
			posmax = posmin;
			posmin = aux;
		}
		//	cout << endl << posmin << " " << posmax;
			//Ponemos huecos en las tareas correspondientes al sector en los hijos inversos
		for (int i = posmin; i <= posmax; i++) {
			int tareaind1 = Poblacion[ind1].vectorSecuenciaOperaciones[i];
			int tareaind2 = Poblacion[ind2].vectorSecuenciaOperaciones[i];
			PoblacionAux[ind1].vectorSecuenciaOperaciones[PoblacionAux[ind1].vectorPosicionOperaciones[tareaind2]] = -1;
			PoblacionAux[ind2].vectorSecuenciaOperaciones[PoblacionAux[ind2].vectorPosicionOperaciones[tareaind1]] = -1;
		}

		//Desplazamos a izquierda
		for (int i = indice1; i < posmin; i++) {
			if (PoblacionAux[ind1].vectorSecuenciaOperaciones[i] == -1) {
				//Si entramos aqui hay que buscar un elemento a la derecha para traerlo aqui
				int traer = i + 1;
				while (PoblacionAux[ind1].vectorSecuenciaOperaciones[traer] == -1) traer++;
				PoblacionAux[ind1].vectorSecuenciaOperaciones[i] = PoblacionAux[ind1].vectorSecuenciaOperaciones[traer];
				PoblacionAux[ind1].vectorSecuenciaEdges[i] = PoblacionAux[ind1].vectorSecuenciaEdges[traer];
				PoblacionAux[ind1].vectorSecuenciaOperaciones[traer] = -1;
			}
			if (PoblacionAux[ind2].vectorSecuenciaOperaciones[i] == -1) {
				//Si entramos aqui hay que buscar un elemento a la derecha para traerlo aqui
				int traer = i + 1;
				while (PoblacionAux[ind2].vectorSecuenciaOperaciones[traer] == -1) traer++;
				PoblacionAux[ind2].vectorSecuenciaOperaciones[i] = PoblacionAux[ind2].vectorSecuenciaOperaciones[traer];
				PoblacionAux[ind2].vectorSecuenciaEdges[i] = PoblacionAux[ind2].vectorSecuenciaEdges[traer];
				PoblacionAux[ind2].vectorSecuenciaOperaciones[traer] = -1;
			}
		}

		//Desplazamos a derecha
		for (int i = limite - 1; i > posmax; i--) {
			if (PoblacionAux[ind1].vectorSecuenciaOperaciones[i] == -1) {
				//Si entramos aqui hay que buscar un elemento a la izquierda para traerlo aqui
				int traer = i - 1;
				while (PoblacionAux[ind1].vectorSecuenciaOperaciones[traer] == -1) traer--;
				PoblacionAux[ind1].vectorSecuenciaOperaciones[i] = PoblacionAux[ind1].vectorSecuenciaOperaciones[traer];
				PoblacionAux[ind1].vectorSecuenciaEdges[i] = PoblacionAux[ind1].vectorSecuenciaEdges[traer];
				PoblacionAux[ind1].vectorSecuenciaOperaciones[traer] = -1;
			}
			if (PoblacionAux[ind2].vectorSecuenciaOperaciones[i] == -1) {
				//Si entramos aqui hay que buscar un elemento a la izquierda para traerlo aqui
				int traer = i - 1;
				while (PoblacionAux[ind2].vectorSecuenciaOperaciones[traer] == -1) traer--;
				PoblacionAux[ind2].vectorSecuenciaOperaciones[i] = PoblacionAux[ind2].vectorSecuenciaOperaciones[traer];
				PoblacionAux[ind2].vectorSecuenciaEdges[i] = PoblacionAux[ind2].vectorSecuenciaEdges[traer];
				PoblacionAux[ind2].vectorSecuenciaOperaciones[traer] = -1;
			}
		}

		//Y Por ultimo rellenamos los sectores en cada hijo
		for (int i = posmin; i <= posmax; i++) {
			PoblacionAux[ind1].vectorSecuenciaOperaciones[i] = Poblacion[ind2].vectorSecuenciaOperaciones[i];
			PoblacionAux[ind2].vectorSecuenciaOperaciones[i] = Poblacion[ind1].vectorSecuenciaOperaciones[i];
			PoblacionAux[ind1].vectorSecuenciaEdges[i] = Poblacion[ind2].vectorSecuenciaEdges[i];
			PoblacionAux[ind2].vectorSecuenciaEdges[i] = Poblacion[ind1].vectorSecuenciaEdges[i];
		}
		//Y para acabar construimos el posicionTareasVector
		for (int i = indice1; i < limite; i++) {
			PoblacionAux[ind1].vectorPosicionOperaciones[PoblacionAux[ind1].vectorSecuenciaOperaciones[i]] = i;
			PoblacionAux[ind2].vectorPosicionOperaciones[PoblacionAux[ind2].vectorSecuenciaOperaciones[i]] = i;
		}

	}

}

//Operador de mutacion
void Mutacion(Individuo &ind, int c, int t) {

	int numM1 = Aleatorio(1, c);

	for (int i = 0; i < numM1; i++) {   // movimientos 1swap
		int numM = Aleatorio(0, 1);
		//Si la mutacion es 1 modificamos la posicion de una puerta ps a otra poscion adyacente
		if (numM == 1) {
			int aleat1 = Aleatorio(itActual*len, itFinal*len - 1);
			//int aleat2 = Aleatorio(0, n*2-1);
			//ind.vectorSecuenciaEdges[aleat1] = aleat2;
			edge ed = aristasDelChip[ind.vectorSecuenciaEdges[aleat1]];
			ind.vectorSecuenciaEdges[aleat1] = edgeAdy[ed.pos1][ed.pos2][Aleatorio(0, edgeAdy[ed.pos1][ed.pos2].size() - 1)];
			//ind.vectorSecuenciaEdges[aleat1] = Aleatorio(0,n*2-1);
		}
		else if (numM == 2) {
			//Si la mutacion es 2 modificamos la posicion de una puerta ps a otra poscion aleatoria
			int aleat1 = Aleatorio(itActual*len, itFinal*len - 1);
			ind.vectorSecuenciaEdges[aleat1] = Aleatorio(0, n * 2 - 1);
		}
		//Si la mutacion es <=0 hacemos ese numero de movimientos swap aleatorios
		else {
			int per = Aleatorio(1, n / 3);
			for (int i = 0; i < per; ++i) 
			{
				int aleat1 = Aleatorio(itActual * len, itFinal * len - 1);
				int aleat2 = Aleatorio((aleat1 / len) * len, (aleat1 / len + 1) * len - 1);
				while (aleat1 == aleat2) aleat2 = Aleatorio((aleat1 / len) * len, (aleat1 / len + 1) * len - 1);
				int auxxxx1 = ind.vectorSecuenciaOperaciones[aleat1];
				int auxxxx2 = ind.vectorSecuenciaEdges[aleat1];
				ind.vectorSecuenciaOperaciones[aleat1] = ind.vectorSecuenciaOperaciones[aleat2];
				ind.vectorSecuenciaEdges[aleat1] = ind.vectorSecuenciaEdges[aleat2];
				ind.vectorPosicionOperaciones[ind.vectorSecuenciaOperaciones[aleat1]] = aleat1;
				ind.vectorSecuenciaOperaciones[aleat2] = auxxxx1;
				ind.vectorSecuenciaEdges[aleat2] = auxxxx2;
				ind.vectorPosicionOperaciones[auxxxx1] = aleat2;
			}
		}
	}
}


bool Iguales(Individuo &ind1, Individuo &ind2) {
	if (ind1.qstatesBaseSuma != ind2.qstatesBaseSuma) return false;
	if (ind1.makespan != ind2.makespan) return false;
	return true;
}

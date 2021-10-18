# include <iostream>
#include <cmath>
# include "QBit.h"
using namespace std;


class ListaQBits{
	private:
		QBit* raiz;
		QBit* verificarAntQStates(ListaQBits* listaQBitsFinal, int nuevoQState);
		QBit* calcularHeuristicas(ListaQBits* qBitsVecinos, ListaQBits* listaQBitsInic, int** distancias);
		float calcularHeuristica(QBit* reco, ListaQBits* qBitsVecinos, int** distancias);
		void imprimirQBit(QBit* qBit);
		void imprimirQBitFuerza(QBit* qBit);
		
		
		
	public:
		ListaQBits();
		void imprimir();
		QBit* getRaiz();
		void crearListaConFuerza(int**);
		void iniciarListaConFuer(int** qBitsOrdInvert, int numeroQbits);
		QBit* extraerPrimero();
		void insertarQBit(QBit* q);
		QBit* extraerHeredero(ListaQBits* listaQBitsInic, ListaQBits* listaQBitsFinal, QState* nuevoQState , operation* puertasDelProblema, int numeroPuertas, int** distancias);	
		void insertarCopQBit(QBit* q);
		QBit* extraerQBit(QBit* qBitAExtr);
		int* devolverListaIni(int tagnoLista);
		void iniciarLista(int* vectorInicios);
		int devolverState(int qBit);
		void imprimirFuerza();
		void InprimirVecinosLog();
		
};

ListaQBits::ListaQBits()
{
    raiz = NULL;
}

int ListaQBits::devolverState(int qBit){
	QBit* reco = raiz;
	while(reco != NULL && reco->getNumero() != qBit){
		reco = reco->getSig();
	}
	return reco->qStateCont->getQState();
}

void ListaQBits::iniciarLista(int* vectorInicios){
	QBit* reco = raiz;
	while(reco != NULL){
		vectorInicios[reco->getNumero()] = reco->qStateCont->getQState();
		reco = reco->getSig();
	}
}

int* ListaQBits::devolverListaIni(int tamagnoLista){
	cout << "tamagnoLista: " << tamagnoLista << endl;
	int* liataADevol = new int[tamagnoLista];
	QBit* reco=raiz;
	while(reco!=NULL){
		cout << "QBit: " <<  reco->getNumero() + 1 << " QState: " << reco->qStateCont->getQState() + 1 << endl;
		liataADevol[reco->getNumero()] = reco->qStateCont->getQState();
		reco = reco->getSig();
	}
	for(int i; i<tamagnoLista; ++i){
		cout << liataADevol[i] << endl;
	}
	return liataADevol;
}

QBit* ListaQBits::extraerHeredero(ListaQBits* listaQBitsInic, ListaQBits* listaQBitsFinal, QState* nuevoQState, operation* puertasDelProblema, int numeroPuertas, int** distancias)
{
	QBit* qBitSig;
	QBit* qBitVec;
	ListaQBits* qBitsVecinos = new ListaQBits();
	qBitSig=NULL;
	for(int i=0; i < numeroPuertas; ++i)
	{
		if(nuevoQState->getQState() == puertasDelProblema[i].q1)
		{
			qBitVec = verificarAntQStates(listaQBitsFinal, puertasDelProblema[i].q2);
			if(qBitVec != NULL)
			{
				qBitsVecinos->insertarCopQBit(qBitVec);
			}
		}
		else if(nuevoQState->getQState() == puertasDelProblema[i].q2)
		{
			qBitVec = verificarAntQStates(listaQBitsFinal, puertasDelProblema[i].q1);
			if(qBitVec != NULL)
			{
				qBitsVecinos->insertarCopQBit(qBitVec);
			}
		}
	}
	if(qBitsVecinos->getRaiz() != NULL)
	{
		cout << "Los vecinos logicos del QState " << nuevoQState->getQState() + 1<< " son:" << endl;
		qBitsVecinos->InprimirVecinosLog();
		cout << endl;
		qBitSig = calcularHeuristicas(qBitsVecinos, listaQBitsInic, distancias);
		qBitSig->qStateCont = nuevoQState;
	}
	return qBitSig;
}

void ListaQBits::InprimirVecinosLog(){
	QBit* reco = raiz;
	if(reco == NULL)
	{
		cout << "No hay vecinos logicos." << endl;
	}
	else
	{
		while(reco != NULL)
		{
			cout << "El QState " << reco->qStateCont->getQState() + 1 << " en el QBit " << reco->getNumero() + 1 << endl;
			reco = reco->getSig();
		}
	}
}


QBit* ListaQBits::verificarAntQStates(ListaQBits* listaQBitsFinal, int nuevoQState){
	QBit* qBitDevolver;
	qBitDevolver = NULL;
	QBit* reco = listaQBitsFinal->getRaiz();
	while(reco != NULL)
	{
		if(reco->qStateCont->getQState() == nuevoQState)
		{
			qBitDevolver = reco;
		}
		reco = reco->getSig();
	}
	return qBitDevolver;
}

QBit* ListaQBits::calcularHeuristicas(ListaQBits* qBitsVecinos, ListaQBits* listaQBitsInic, int** distancias)
{
	float mejorHeur = 0;
	float heur = 0;
	QBit* qBitSig;
	QBit* reco;
	qBitSig = NULL;
	reco = listaQBitsInic->getRaiz();
	cout << "Calculamos las heuristicas para los QBits libres:" << endl;
	while(reco != NULL)
	{
		heur = calcularHeuristica(reco,qBitsVecinos, distancias);
		if(heur > mejorHeur)
		{
			mejorHeur=heur;
			qBitSig = reco;
			cout << "Actualizamos QBit." << endl;
		}
		reco = reco->getSig();
	}
	return qBitSig;
}

float ListaQBits::calcularHeuristica(QBit* nuevoQBit, ListaQBits* qBitsVecinos, int** distancias)
{
	double valorHeur=0;
	int sumatorio = 0;
	QBit* reco = qBitsVecinos->getRaiz();
	while(reco != NULL)
	{
		sumatorio = sumatorio + distancias[nuevoQBit->getNumero()][reco->getNumero()];
		reco = reco->getSig();
	}
	
	valorHeur = (float)nuevoQBit->getFuerzaConex()/(float)sumatorio;
	cout << "La heuristica del QBit " << nuevoQBit->getNumero() + 1 << " es " << valorHeur << endl;
	return valorHeur;
}

QBit* ListaQBits::getRaiz()
{
	return raiz;
}

void ListaQBits::iniciarListaConFuer(int** qBitsOrdInvert, int numeroQbits)
{
	QBit* reco;
	QBit* nuevo;
	for(int i = 0; i < numeroQbits; ++i)
	{
		reco = raiz;
		nuevo = new QBit(qBitsOrdInvert[i][0],reco);
		nuevo->setFuerzaConexion(qBitsOrdInvert[i][1]);
		raiz = nuevo;
	}
}

QBit* ListaQBits::extraerPrimero()
{
	QBit* extraer = raiz;
	raiz = extraer->getSig();
	return extraer;
}

QBit* ListaQBits::extraerQBit(QBit* qBitAExtr){
	QBit* reco;
	QBit* ant;
	reco = raiz;
	if(reco == NULL){
		cout << "No es posible extraer QBit, ya que la lista esta vacia. " << endl;
	}else if(reco == qBitAExtr){
		raiz = raiz->getSig();
	}
	else
	{
		ant = reco;
		reco = reco->getSig();
		while(reco != NULL && reco != qBitAExtr){
			ant = reco;
			reco = reco->getSig();
		}
		if(reco == NULL){
			cout << "No es posible extraer QBit, ya que este no esta en la lista. " << endl;
		}else{
			ant->setSig(reco->getSig());
		}
	}
	return reco;
}

void ListaQBits::insertarQBit(QBit* q)
{
	if(raiz == NULL){
		q->setSig(raiz);
		raiz = q;
	}else{
		QBit* reco = raiz;
		while(reco->getSig() != NULL)
		{
			reco = reco->getSig();
		}
		q->setSig(reco->getSig());
		reco->setSig(q);
	}
}

void ListaQBits::insertarCopQBit(QBit* q)
{
	QBit* copQBit = new QBit(q->getNumero());
	copQBit->qStateCont = q->qStateCont;
	if(raiz == NULL){
		copQBit->setSig(raiz);
		raiz = copQBit;
	}else{
		QBit* reco = raiz;
		while(reco->getSig() != NULL)
		{
			reco = reco->getSig();
		}
		copQBit->setSig(reco->getSig());
		reco->setSig(copQBit);
	}
}

void ListaQBits::imprimir()
{
    QBit* reco = raiz;
    cout << "Resumen Inicializaciones:" << endl;
    while (reco != NULL)
    {
        imprimirQBit(reco);
        reco = reco->getSig();
    }
    cout << "\n";
}

void ListaQBits::imprimirQBit(QBit* qBit)
{
	cout << "Inicializamos el QState ";
	if(qBit->qStateCont != NULL){
		cout << qBit->qStateCont->getQState() + 1;
	}else{
		cout << "nulo";
	}
	cout << " en el QBit " << qBit->getNumero() + 1 << "\n";
}

void ListaQBits::imprimirFuerza()
{
	cout << "Lista de QBits segun su fuerza de conexion: " << endl;
    QBit* reco = raiz;
    int cont = 1;
    while (reco != NULL)
    {
    	cout << cont++ << " : ";
        imprimirQBitFuerza(reco);
        reco = reco->getSig();
    }
    cout << "\n";
}

void ListaQBits::imprimirQBitFuerza(QBit* qBit)
{
	cout << "El QBit ";
	cout << qBit->getNumero();
	cout << " tiene fuerza ";
	cout << qBit->getFuerzaConex() << endl;
}









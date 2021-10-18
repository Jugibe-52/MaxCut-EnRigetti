# include <iostream>
#include "../funciones/calcularOrdenQstates.h"
# include "QState.h"
using namespace std;

class ListaQState
{
private:
    QState* raiz;
public:
    ListaQState();
    void iniciarLista(int* lista, int numeroQStates);
    void insertarQState(int qState);
    void insertarQState(QState* qState);
    void insertarLista(ListaQState* lista);
    void imprimir();
    class QState* getRaiz();
    QState* extraerPrimero(); 
};

ListaQState::ListaQState()
{
    raiz = NULL;
}

class QState* ListaQState::getRaiz()
{
	return raiz;
}

void ListaQState::iniciarLista(int* lista, int numeroQStates)
{
	QState* nuevo;
	QState* reco;
	
	for(int i = 0; i < numeroQStates; ++i)
	{
		reco = raiz;
		nuevo = new QState(lista[i],reco);
		raiz = nuevo;
	}
}

void ListaQState::insertarQState(int x)
{
    QState* nuevo = new QState();
    nuevo->setQState(x);
    if (raiz == NULL) 
    {
        raiz = nuevo;
    }
    else 
    {
    	QState *ant = raiz;
        QState *reco = ant->getSig();
        
        while (reco != NULL) 
        {
        	ant->setSig(reco);
        	reco->setSig(reco->getSig());
        }
        ant->setSig(nuevo);
        nuevo->setSig(reco);
    }
}

void ListaQState::insertarQState(QState* q)
{
	if(raiz = NULL){
		q->setSig(raiz);
		raiz = q;
	}else{
		QState* reco = raiz;
		while(reco->getSig() != NULL)
		{
			reco = reco->getSig();
		}
		q->setSig(reco->getSig());
		reco->setSig(q);
	}
}

void ListaQState::insertarLista(ListaQState* lista)
{
	QState *cabeza = raiz;
	while(cabeza != NULL)
	{
		cabeza = cabeza->getSig();
	}
	cabeza->setSig(lista->getRaiz());
}

void ListaQState::imprimir()
{
    QState *reco = raiz;
    int contador = 0;
    cout << "Orden Seleccion QStates:\n";
    while (reco != NULL)
    {
        cout << contador++ << ": "<< reco->getQState() + 1 << endl;
        reco = reco->getSig();
    }
    cout << "\n";
}

QState* ListaQState::extraerPrimero(){
	QState* extraer = raiz;
	raiz = extraer->getSig();
	return extraer;
}












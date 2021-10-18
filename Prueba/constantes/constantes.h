#ifndef PCH_H
#define PCH_H
#endif //PCH_H

// int ITERATIONS;
#define INICIALES 0 //Inicializacion de qstates en distintas posiciones

#define DEBUG 1
#define DEBUG2 0

#define MUYNEGATIVO -10000000
#define MUYPOSITIVO 10000000
#define MAXVECINOS 100000
#define MAXIMOMAKESPAN 1000000

#define METODO 1    // 0=Scatter Search, 1=Genetico Generacional+BL, 2=Tabu Search sola


#define TASGETIREN 1

//Para controlar el maximo gradiente
#define VECINDADMAXIMOGRADIENTE 4

//Para controlar el scatter search
#define MAXITERSINMEJORARSCATTERSEARCH 100      // <------ esto se usa en scatter search (iteraciones)
#define TAMANOPOBLACION 10      //Esto se usa para scatter search
#define NUMMEJORES 10             //Esto se usa en RefSetUpdate, para determinar el tamano del subconjunto de mejores y del subconjunto de diversos
#define TAMANOSOLUCIONESINICIALES 20     //Esto se usa solo para scatter search
#define NUMMEJORESINICIAL 5       //Esto define, para generar Poblacion, cuantos de los mejores se eligen del conjunto P, el resto se rellenara con diversos
#define UTILIZARDIVERSIFICACION 1        //Si esta a -1 entonces en cuanto vaya a hacer una diversificacion se acaba
#define CUANTOSMEJORESDIVERSIFICACION 1     //Esto se usa en CrearPoblacionInicialDiversificacion, para ver con cuantos de los mejores de la poblacion nos quedamos
#define DISTANCIAMINIMAENTREELEMENTOS 20

//Para controlar el genetico
#define SOLOBL 1
// int TAMANOPOBLACIONGENETICO;
#define MAXGENERACIONESSINMEJORAGENETICO 800 //300//800 //Condicion de parada, numero de generaciones sin mejora 


#define PRIMEROS 10000 //Numero de individuos que se generan de forma heuritica al crear la poblacion
#define SELECCION 3 // 1=los dos mejores, 2=mejor y peor, 3=dos mejores con diferente makespan, 4=dos mejores que sean diferentes
#define TIPOCRUCE 0 // 0=JOX  1=NWOX
#define PROBCRUCE 100  //Esto sera un numero entre 0 y 100
#define ETAPAS 1 //Si esta a 1 realiza la busqueda por etapas, tratando solo una iteracion de cada vez
#define ETAPASF 0 //Si esta a 1 y realiza la busqueda por etapas, hace una busqueda final con todo
#define TIPOMUTACION 1 // numero de mutaciones maximo en la mutacion
#define TIPOMUTACIOND 5 // numero de mutaciones maximo en la diversificacion
#define PROBMUTACION 40 //Esto sera un numero entre 0 y 100
#define MUTACIONADAPTATIVA 0  // Aumentara PORCENTMUTACIONADAPTATIVA% si hay una generacion sin mejoras. Volvera a 0% si hay una con mejoras.
#define PORCENTMUTACIONADAPTATIVA 1 
#define ELIMINARIGUALES 10   //Cada cuantas generaciones sin mejora utiliza el operador de diversificacion
#define ELIMINARIGUALESPOR 60 //Porcentaje en que se realiza la diversificacion
#define NUMISLANDS 1


#define EPOCH 20   // Esto indica cada cuantas generaciones se intercambian individuos entre los distintos islands

//Para solo busqueda tabu unica y a tiempo, pon solobusquedatabu a 1, tamanosolucionesiniciales 1, maxiteracionestotales y sin mejora a -1, y ejecuta el programa como pruebas segundos 100000000 azar

//Usados por la busqueda tabu
#define TAMANOSOLUCIONESINICIALESTS 40   //Utilizado solo si se ejecuta la TS por separado
#define VECINDADBUSQUEDATABUALEATORIA 0
#define VECINDADBUSQUEDATABU 4                 //vecindad 0 Swap vecindad 1 Insert vecindad 2 2BlockInsert 3 3BlockInsert 4 BlockMove 5 BlockSwap 6 BlockMove+BlockSwap 7 BlockReverse 8 BlockMove+BlockSwap+BlockReverse
#define TAMANO_MIN_INICIAL 4
#define TAMANO_MAX_INICIAL 12
#define LONGITUD_INICIAL_LISTA_TABU 8
#define MAX_ITERACIONES_TOTALES -1
#define MAX_ITERACIONES_SIN_MEJORA 50
#define RESTABLECER_MIN_MAX true
#define NUM_ITERACIONES_RESTABLECER_MIN_MAX 20
//#define TAMANO_DINAMICO true
#define MAX_MOVIMIENTOS_CICLICOS 4             // -1 para no utilizar deteccion de ciclos


//Usados por el path relinking
#define VECINDADPATHRELINKING 0
#define ELEGIRMEDIOTRAYECTORIA 0
#define GUARDARTRAYECTORIAPATHRELINKING 1
#define ELEGIRVARIOSTRAYECTORIA 1     // Para esto debe estar guardartrayectoria a 1 y elegirmediotrayectoria a 0
#define CUANTOSVARIOSTRAYECTORIA 3


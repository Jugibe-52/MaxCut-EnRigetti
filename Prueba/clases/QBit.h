#include <iostream>
using namespace std;

class QBit
{
	private:
		int numero;
		int fuerzaConex;
		
		QBit* sig;
	public:
		
		QBit(){
		}
		QBit(int n){
			numero = n;
			sig = NULL;
			qStateCont = NULL;
		}
		QBit(int n,QBit* s){
			numero = n;
			sig = s;
			qStateCont = NULL;
		}
		
		QState* qStateCont;
		
		int getNumero(){
			return numero;
		}
		int getFuerzaConex(){
			return fuerzaConex;
		}
		QBit* getSig(){
			return sig;
		}
		void setNumero(int n){
			numero = n;
		}
		void setSig(QBit* s){
			sig = s;
		}
		void setFuerzaConexion(int f){
			fuerzaConex = f;
		}
};

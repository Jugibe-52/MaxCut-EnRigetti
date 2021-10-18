#include <iostream>
using namespace std;

class QState
{
	private:
		int qState;
		QState* sig;
	public:
		QState(){
		}
		QState(int q)
		{
			qState = q;
			sig = NULL;
		}
		QState(int q, QState* s)
		{
			qState = q;
			sig = s;
		}
		void setQState(int q)
		{
			qState = q;
		}
		void setSig(QState *s)
		{
			sig = s;
		}
		int getQState()
		{
			return qState;
		}
		class QState* getSig()
		{
			return sig;
		}
};


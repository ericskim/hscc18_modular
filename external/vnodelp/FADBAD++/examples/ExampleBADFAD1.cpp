#include "badiff.h"
#include "fadiff.h"
#include <iostream>
using namespace std;

B< F<double> > func(const B< F<double> >& x, const B< F<double> >& y)
{
	B< F<double> > z=sqrt(x);
	return y*z+sin(z);
}

int main()
{
	B< F<double> > x,y,f;      // Declare variables x,y,f
	x=1;                       // Initialize variable x
	y=2;                       // Initialize variable y
	x.x().diff(0,2);           // Second order wrt. x
	y.x().diff(1,2);           // Second order wrt. y
	f=func(x,y);               // Evaluate function and record DAG
	f.diff(0,1);               // Differentiate f
	double fval=f.x().x();     // Value of function
	double dfdx=x.d(0).x();    // Value of df/dx
	double dfdy=y.d(0).x();    // Value of df/dy
	double dfdxdx=x.d(0).d(0); // Value of df²/dx²
	double dfdxdy=x.d(0).d(1); // Value of df²/dxdy
	double dfdydx=y.d(0).d(0); // Value of df²/dydx
	double dfdydy=y.d(0).d(1); // Value of df²/dy²

	cout << "f(x,y)=" << fval << endl;
	cout << "df/dx(x,y)=" << dfdx << endl;
	cout << "df/dy(x,y)=" << dfdy << endl;
	cout << "df²/dx²(x,y)=" << dfdxdx << endl;
	cout << "df²/dxdy(x,y)=" << dfdxdy << endl;
	cout << "df²/dydx(x,y)=" << dfdydx << endl;
	cout << "df²/dy²(x,y)=" << dfdydy << endl;

	return 0;
}

/*
 * RungeKutta4.hh
 *
 *  created on: 22.04.2015
 *      author: rungger
 */

#ifndef RUNGEKUTTA4_HH_
#define RUNGEKUTTA4_HH_

/* class: RungeKutta4 
 * Fixed step size ode solver implementing a Runge Kutta scheme of order 4 */
class OdeSolver {
public:
  /* function: ()
   *
   * Input:
   * rhs - rhs of ode  dx/dt = rhs(x,u)
   * x - current state x
   * u - current input u
   * dim - state space dimension
   * tau - sumpling time
   * nint - number of intermediate steps (default = 10)
   */
  template<class RHS, class X, class U>
  void operator()(RHS rhs, X &x, U &u, const int dim, const double tau, const int nint=10) {
		X k[4];
		X tmp;

    double h=tau/(double)nint;

		for(int t=0; t<nint; t++) {
			rhs(k[0],x,u);
			for(int i=0;i<dim;i++)
				tmp[i]=x[i]+h/2*k[0][i];

			rhs(k[1],tmp, u);
			for(int i=0;i<dim;i++)
				tmp[i]=x[i]+h/2*k[1][i];

			rhs(k[2],tmp, u);
			for(int i=0;i<dim;i++)
				tmp[i]=x[i]+h*k[2][i];

			rhs(k[3],tmp, u);
			for(int i=0; i<dim; i++)
				x[i] = x[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
		}
	}
};

#endif /* RUNGEKUTTA4_HH_ */

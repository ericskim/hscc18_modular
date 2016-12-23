#include "badiff.h" // Backwards automatic differentiation
#include "fadiff.h" // Forwards automatic differentiation
#include "ffadiff.h" // Fast-Forward automatic differentiation

#include "extra/ndf.h"
#include "extra/ndfad.h"

#include <vector>
#include "TestFADBAD.h"

// Functions to differentiate:

struct Add
{
	static const char* Name(){ return "simple addition"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]+v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
struct Add1
{
	static const char* Name(){ return "simple addition 1"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]+117;
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Add2
{
	static const char* Name(){ return "simple addition 2"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=118+v[0];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Sub
{
	static const char* Name(){ return "simple subtraction"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]-v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
struct Sub1
{
	static const char* Name(){ return "simple subtraction 1"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]-117;
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Sub2
{
	static const char* Name(){ return "simple subtraction 2"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=118-v[0];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Mul
{
	static const char* Name(){ return "simple multiplication"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]*v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
struct Mul1
{
	static const char* Name(){ return "simple multiplication 1"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]*117;
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Mul2
{
	static const char* Name(){ return "simple multiplication 2"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=118*v[0];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Div
{
	static const char* Name(){ return "simple division"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]/v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
struct Div1
{
	static const char* Name(){ return "simple division 1"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=v[0]/117;
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Div2
{
	static const char* Name(){ return "simple division 2"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=118/v[0];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Pow
{
	static const char* Name(){ return "simple power"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=pow(v[0],3);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Pow2
{
	static const char* Name(){ return "power"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=pow(v[0],v[1]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[0]=3.3;
		return in;
	}
};
struct Sqr
{
	static const char* Name(){ return "simple square"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=sqr(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Exp
{
	static const char* Name(){ return "simple exponential"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=exp(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Log
{
	static const char* Name(){ return "simple logarithmic"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=log(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Sqrt
{
	static const char* Name(){ return "simple square-root"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=sqrt(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Sin
{
	static const char* Name(){ return "simple sine"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=sin(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Cos
{
	static const char* Name(){ return "simple cosine"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=cos(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct Tan
{
	static const char* Name(){ return "simple tangent"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=tan(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		return in;
	}
};
struct ASin
{
	static const char* Name(){ return "simple arc sine"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=asin(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.5;
		return in;
	}
};
struct ACos
{
	static const char* Name(){ return "simple arc cosine"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=acos(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.5;
		return in;
	}
};
struct ATan
{
	static const char* Name(){ return "simple arc tangent"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=atan(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.5;
		return in;
	}
};
struct Nd
{
	static const char* Name(){ return "normal distribution function"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=nd(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.5;
		return in;
	}
};
struct Cnd
{
	static const char* Name(){ return "cumulative normal distribution function"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=cnd(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.5;
		return in;
	}
};
struct Icnd
{
	static const char* Name(){ return "inverse cumulative normal distribution function"; }
	static const int m=1;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=icnd(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=0.8;
		return in;
	}
};



struct csmap // The cos-sin map.
{
	static const char* Name(){ return "Cos-Sine map"; }
	static const int m=2;
	static const int n=2;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=cos(v[0]+4*v[1]);
		out[1]=sin(4*v[0]+v[1]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
struct func1
{
	static const char* Name(){ return "atan[(a sin x)/(1-a cos x)]"; }
	static const int m=2;
	static const int n=1;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		T x(v[0]);
		T a(v[1]);
		out[0]=atan((a*sin(x))/(1-a*cos(x)));
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(m);
		in[0]=1.5;
		in[1]=-.456;
		return in;
	}
};



// The classes Fdiff<C> and Bdiff are capable of differentiating a function
// f = [f0,f1,..,f(m-1)] : T^n->T^m by using either the forward or the backward
// modes of automatic differentiation. Where T is real, interval, etc.
// The argument when using the evaluating operation on Fdiff is an
// n-dimensional std::vector x in which the function should be evaluated
// and differentiated. The return value is a std::vector that should
// be interpreted as follows:
//		Vout[0..m-1] is the function value.
//		Vout[m..m+n-1] is the value of [df1/dx1, ..., df1/dxn]
//		Vout[m+i*n..m+(i+1)*n-1] is the value of [dfi/dx1, ..., dfi/dxn].

template <class C>
struct Fdiff
{
	static const int m=C::m;
	static const int n=C::n+C::m*C::n;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j,inSize=v.size();;
		std::vector< F<T> > out;
		std::vector< F<T> > in(inSize);
		for(i=0;i<inSize;++i)
		{
			in[i]=v[i];
			in[i].diff(i,inSize);
		}
		out=C()(in);
		int outSize=out.size();
		std::vector<T> retval(outSize*(1+inSize));
		for(j=0;j<outSize;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<inSize;++i)
			{
				retval[outSize+j*inSize+i]=out[j].d(i);
			}
		}
		return retval;
	}
};

template <class C>
struct FFdiff // fast-forward variant; the number of independent vars. is known at compile-time
{
	static const int m=C::m;
	static const int n=C::n+C::m*C::n;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j;
		std::vector< FF<T,C::m> > out;
		std::vector< FF<T,C::m> > in(C::m);
		for(i=0;i<C::m;++i)
		{
			in[i]=v[i];
			in[i].diff(i);
		}
		out=C()(in);
		std::vector<T> retval(C::n*(1+C::m));
		for(j=0;j<C::n;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<C::m;++i)
			{
				retval[C::n+j*C::m+i]=out[j].d(i);
			}
		}
		return retval;
	}
};

template <class C>
struct Bdiff
{
	static const int m=C::m;
	static const int n=C::n+C::m*C::n;
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j,inSize=v.size();
		std::vector< B<T> > out;
		std::vector< B<T> > in(inSize);
		for(i=0;i<inSize;++i)
		{
			in[i]=v[i];
		}
		out=C()(in);
		int outSize=out.size();
		for(j=0;j<outSize;++j)
		{
			out[j].diff(j,outSize);
		}
		std::vector<T> retval(outSize*(1+inSize));
		for(j=0;j<outSize;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<inSize;++i)
			{
				retval[outSize+j*inSize+i]=in[i].d(j);
			}
		}
		return retval;
	}
};

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
	os << "[";
	for(int i=0;i<v.size();i++)
	{
		os << v[i];
		if (i<v.size()-1) os << ",";
	}
	os << "]";
	return os;
}

template <class T>
struct check
{
	void runFB(IReportLog& log);
	void runFFB(IReportLog& log);
};

bool diff(const std::vector<double> v1, const std::vector<double> v2)
{
//	std::cout<<"COMPARE"<<std::endl;
//	std::cout<<v1<<std::endl;
//	std::cout<<v2<<std::endl;

	if (v1.size()!=v2.size()) return true;
	for(unsigned int i=0;i<v1.size();++i)
		if (fabs(v1[i]-v2[i])>1e-6) return true;
	return false;
}

template <class T>
void check<T>::runFB(IReportLog& log)
{
	std::vector<double> v1;
	std::vector<double> v2;
	std::vector<double> v3;

	{
		Fdiff< T > Fd;
		v1=Fd(T::point());
	}

	{
		Bdiff< T > Bd;
		v2=Bd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 1. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 1. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff< Fdiff < T > > FFd;
		v1=FFd(T::point());
	}

	{
		Bdiff< Fdiff < T > > BFd;
		v2=BFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff< Bdiff < T > > FBd;
		v2=FBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff< Bdiff < T > > BBd;
		v2=BBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff< Fdiff < T > > > FFFd;
		v1=FFFd(T::point());
	}

	{
		Bdiff < Fdiff< Fdiff < T > > > BFFd;
		v2=BFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff< Fdiff < T > > > FBFd;
		v2=FBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< Fdiff < T > > > BBFd;
		v2=BBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff< Bdiff < T > > > FFBd;
		v2=FFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff< Bdiff < T > > > BFBd;
		v2=BFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff< Bdiff < T > > > FBBd;
		v2=FBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< Bdiff < T > > > BBBd;
		v2=BBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;
/* Uncomment this if you want to wait even longer for compilation to finish
	{
		Fdiff < Fdiff < Fdiff< Fdiff < T > > > > FFFFd;
		v1=FFFFd(T::point());
	}

	{
		Bdiff < Fdiff < Fdiff< Fdiff < T > > > > BFFFd;
		v2=BFFFd(T::point());
	}
	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Fdiff< Fdiff < T > > > > FBFFd;
		v2=FBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Fdiff< Fdiff < T > > > > BBFFd;
		v2=BBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Bdiff< Fdiff < T > > > > FFBFd;
		v2=FFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Bdiff< Fdiff < T > > > > BFBFd;
		v2=BFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Bdiff< Fdiff < T > > > > FBBFd;
		v2=FBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< Fdiff < T > > > > BBBFd;
		v2=BBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Fdiff< Bdiff < T > > > > FFFBd;
		v2=FFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Fdiff< Bdiff < T > > > > BFFBd;
		v2=BFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Fdiff< Bdiff < T > > > > FBFBd;
		v2=FBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Fdiff< Bdiff < T > > > > BBFBd;
		v2=BBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Bdiff< Bdiff < T > > > > FFBBd;
		v2=FFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Bdiff< Bdiff < T > > > > BFBBd;
		v2=BFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Bdiff< Bdiff < T > > > > FBBBd;
		v2=FBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< Bdiff < T > > > > BBBBd;
		v2=BBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;
*/
}

template <class T>
void check<T>::runFFB(IReportLog& log)
{
	std::vector<double> v1;
	std::vector<double> v2;
	std::vector<double> v3;

	{
		FFdiff< T > Fd;
		v1=Fd(T::point());
	}

	{
		Bdiff< T > Bd;
		v2=Bd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 1. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 1. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff< FFdiff < T > > FFd;
		v1=FFd(T::point());
	}

	{
		Bdiff< FFdiff < T > > BFd;
		v2=BFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff< Bdiff < T > > FBd;
		v2=FBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff< Bdiff < T > > BBd;
		v2=BBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < FFdiff< FFdiff < T > > > FFFd;
		v1=FFFd(T::point());
	}

	{
		Bdiff < FFdiff< FFdiff < T > > > BFFd;
		v2=BFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff< FFdiff < T > > > FBFd;
		v2=FBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< FFdiff < T > > > BBFd;
		v2=BBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < FFdiff< Bdiff < T > > > FFBd;
		v2=FFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < FFdiff< Bdiff < T > > > BFBd;
		v2=BFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff< Bdiff < T > > > FBBd;
		v2=FBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< Bdiff < T > > > BBBd;
		v2=BBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;
/* Uncomment this if you want to wait even longer for compilation to finish
	{
		FFdiff < FFdiff < FFdiff< FFdiff < T > > > > FFFFd;
		v1=FFFFd(T::point());
	}

	{
		Bdiff < FFdiff < FFdiff< FFdiff < T > > > > BFFFd;
		v2=BFFFd(T::point());
	}
	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff < FFdiff< FFdiff < T > > > > FBFFd;
		v2=FBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < FFdiff< FFdiff < T > > > > BBFFd;
		v2=BBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < FFdiff < Bdiff< FFdiff < T > > > > FFBFd;
		v2=FFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < FFdiff < Bdiff< FFdiff < T > > > > BFBFd;
		v2=BFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff < Bdiff< FFdiff < T > > > > FBBFd;
		v2=FBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< FFdiff < T > > > > BBBFd;
		v2=BBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < FFdiff < FFdiff< Bdiff < T > > > > FFFBd;
		v2=FFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < FFdiff < FFdiff< Bdiff < T > > > > BFFBd;
		v2=BFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff < FFdiff< Bdiff < T > > > > FBFBd;
		v2=FBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < FFdiff< Bdiff < T > > > > BBFBd;
		v2=BBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < FFdiff < Bdiff< Bdiff < T > > > > FFBBd;
		v2=FFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < FFdiff < Bdiff< Bdiff < T > > > > BFBBd;
		v2=BFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		FFdiff < Bdiff < Bdiff< Bdiff < T > > > > FBBBd;
		v2=FBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< Bdiff < T > > > > BBBBd;
		v2=BBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;
*/
}



void TestFADBAD::run(IReportLog& rlog)
{
	{ check<Add> chk; chk.runFB(rlog); }
	{ check<Add1> chk; chk.runFB(rlog); }
	{ check<Add2> chk; chk.runFB(rlog); }
	{ check<Sub> chk; chk.runFB(rlog); }
	{ check<Sub1> chk; chk.runFB(rlog); }
	{ check<Sub2> chk; chk.runFB(rlog); }
	{ check<Mul> chk; chk.runFB(rlog); }
	{ check<Mul1> chk; chk.runFB(rlog); }
	{ check<Mul2> chk; chk.runFB(rlog); }
	{ check<Div> chk; chk.runFB(rlog); }
	{ check<Div1> chk; chk.runFB(rlog); }
	{ check<Div2> chk; chk.runFB(rlog); }
	{ check<Pow> chk; chk.runFB(rlog); }
	{ check<Pow2> chk; chk.runFB(rlog); }
	{ check<Sqr> chk; chk.runFB(rlog); }
	{ check<Exp> chk; chk.runFB(rlog); }
	{ check<Log> chk; chk.runFB(rlog); }
	{ check<Sqrt> chk; chk.runFB(rlog); }
	{ check<Sin> chk; chk.runFB(rlog); }
	{ check<Cos> chk; chk.runFB(rlog); }
	{ check<Tan> chk; chk.runFB(rlog); }
	{ check<ASin> chk; chk.runFB(rlog); }
	{ check<ACos> chk; chk.runFB(rlog); }
	{ check<ATan> chk; chk.runFB(rlog); }
	{ check<Nd> chk; chk.runFB(rlog); }
	{ check<Cnd> chk; chk.runFB(rlog); }
	{ check<Icnd> chk; chk.runFB(rlog); }
	{ check<csmap> chk; chk.runFB(rlog); }
	{ check<func1> chk; chk.runFB(rlog); }

	{ check<Add> chk; chk.runFFB(rlog); }
	{ check<Add1> chk; chk.runFFB(rlog); }
	{ check<Add2> chk; chk.runFFB(rlog); }
	{ check<Sub> chk; chk.runFFB(rlog); }
	{ check<Sub1> chk; chk.runFFB(rlog); }
	{ check<Sub2> chk; chk.runFFB(rlog); }
	{ check<Mul> chk; chk.runFFB(rlog); }
	{ check<Mul1> chk; chk.runFFB(rlog); }
	{ check<Mul2> chk; chk.runFFB(rlog); }
	{ check<Div> chk; chk.runFFB(rlog); }
	{ check<Div1> chk; chk.runFFB(rlog); }
	{ check<Div2> chk; chk.runFFB(rlog); }
	{ check<Pow> chk; chk.runFFB(rlog); }
	{ check<Pow2> chk; chk.runFFB(rlog); }
	{ check<Sqr> chk; chk.runFFB(rlog); }
	{ check<Exp> chk; chk.runFFB(rlog); }
	{ check<Log> chk; chk.runFFB(rlog); }
	{ check<Sqrt> chk; chk.runFFB(rlog); }
	{ check<Sin> chk; chk.runFFB(rlog); }
	{ check<Cos> chk; chk.runFFB(rlog); }
	{ check<Tan> chk; chk.runFFB(rlog); }
	{ check<ASin> chk; chk.runFFB(rlog); }
	{ check<ACos> chk; chk.runFFB(rlog); }
	{ check<ATan> chk; chk.runFFB(rlog); }
	{ check<Nd> chk; chk.runFB(rlog); }
	{ check<Cnd> chk; chk.runFB(rlog); }
	{ check<Icnd> chk; chk.runFB(rlog); }
	{ check<csmap> chk; chk.runFFB(rlog); }
	{ check<func1> chk; chk.runFFB(rlog); }
}

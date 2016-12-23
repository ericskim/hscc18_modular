/*127:*/
#line 83 "./matrix.w"

#ifndef VECTOR_MATRIX
#define VECTOR_MATRIX

#include <vector> 
#include "vnodeinterval.h"
using namespace std;
using namespace v_bias;
namespace v_blas{
/*122:*/
#line 11 "./matrix.w"


/*80:*/
#line 118 "./integfun.w"


#include <vector> 
#include "vnodeinterval.h"
using namespace std;
using namespace v_bias;
typedef vector<interval> iVector;






/*:80*/
#line 13 "./matrix.w"


typedef vector<double> pVector;
typedef vector<vector<double> > pMatrix;
typedef vector<vector<interval> > iMatrix;




/*:122*/
#line 92 "./matrix.w"

/*124:*/
#line 44 "./matrix.w"


template<class Matrix> 
inline void sizeM(Matrix&A,unsigned int n)
{
A.resize(n);
for(unsigned int i= 0;i<A.size();i++)
A[i].resize(n);
}


template<class Matrix> 
inline unsigned int sizeM(const Matrix&A)
{
return A.size();
}

template<class Vector> 
inline void sizeV(Vector&a,unsigned int n)
{
a.resize(n);
}




template<class Vector> inline unsigned int sizeV(const Vector&a)
{
return a.size();
}



/*:124*/
#line 93 "./matrix.w"

}
#endif



#line 1 "./linalg.w"

/*:127*/

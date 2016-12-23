/*291:*/
#line 1608 "./iho.w"

#include <cstdlib> 
#include <cassert> 
#include "vnodeinterval.h"
#include "basiclinalg.h"

namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*285:*/
#line 24 "./sortcols.w"

struct index_norm{int index;double norm;};


/*:285*/
#line 1617 "./iho.w"

/*288:*/
#line 61 "./sortcols.w"


inline int v_compare(const void*a1,const void*b1)
{
index_norm*a= (index_norm*)a1;
index_norm*b= (index_norm*)b1;

if(a->norm> b->norm)
return-1;
if(a->norm<b->norm)
return 1;

return 0;
}


/*:288*/
#line 1618 "./iho.w"

/*284:*/
#line 8 "./sortcols.w"


void sortColumns(pMatrix&B,const
pMatrix&A)

{
/*286:*/
#line 31 "./sortcols.w"

int n= v_blas::sizeM(B);
index_norm*b= new index_norm[n];
pVector tmp;

v_blas::sizeV(tmp,n);
for(int j= 0;j<n;j++)
{
b[j].index= j;
getColumn(tmp,A,j);
b[j].norm= v_blas::norm2(tmp);
}

/*:286*/
#line 14 "./sortcols.w"


/*287:*/
#line 45 "./sortcols.w"

bool sorting_needed= false;
for(int i= 0;i<n-1;i++)
{
if(b[i].norm<b[i+1].norm)
{
sorting_needed= true;
break;
}
}



/*:287*/
#line 16 "./sortcols.w"


/*289:*/
#line 80 "./sortcols.w"

if(sorting_needed)
std::qsort((void*)b,n,sizeof(index_norm),v_compare);

for(int j= 0;j<n;j++)
{
getColumn(tmp,A,b[j].index);
setColumn(B,tmp,j);
}
delete[]b;








#line 1602 "./iho.w"



/*:289*/
#line 18 "./sortcols.w"

}

/*:284*/
#line 1619 "./iho.w"

}



/*:291*/

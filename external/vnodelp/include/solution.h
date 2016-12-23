/*184:*/
#line 157 "./solut.w"

#ifndef SOLUTION_H
#define SOLUTION_H
namespace vnodelp{
using namespace v_bias;
using namespace v_blas;

/*178:*/
#line 44 "./solut.w"

class Solution
{
public:
Solution(int n);
void init(const v_bias::interval&t0,const iVector&y0);
v_bias::interval t;
pVector u;
iVector y;
iVector alpha,r,rQR;
pMatrix S,A,Q;
};


/*:178*/
#line 164 "./solut.w"

/*181:*/
#line 132 "./solut.w"

class Apriori
{
public:
Apriori(int n);
void init(const v_bias::interval t0,const iVector y0);
v_bias::interval t;
iVector y;
};



/*:181*/
#line 165 "./solut.w"

}
#endif

/*:184*/

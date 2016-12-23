/*185:*/
#line 170 "./solut.w"

#include "vnodeinterval.h"
#include "basiclinalg.h"
#include "solution.h"
#include "intvfuncs.h"


namespace vnodelp{
/*179:*/
#line 60 "./solut.w"

Solution::Solution(int n)
{
sizeV(u,n);
sizeV(alpha,n);
sizeV(r,n);
sizeV(rQR,n);
sizeV(y,n);

sizeM(S,n);
sizeM(A,n);
sizeM(Q,n);
}


/*:179*/
#line 178 "./solut.w"

/*180:*/
#line 91 "./solut.w"

void Solution::init(const v_bias::interval&t0,const iVector&y0)
{

t= t0;
y= y0;
midpoint(u,y0);

subViVp(alpha,y0,u);


setV(r,0.0);
setV(rQR,0.0);

setId(S);
setId(Q);
setId(A);


}





/*:180*/
#line 179 "./solut.w"

/*182:*/
#line 145 "./solut.w"

Apriori::Apriori(int n)
{sizeV(y,n);}

/*:182*/
#line 180 "./solut.w"

/*183:*/
#line 150 "./solut.w"

void Apriori::init(const v_bias::interval t0,const iVector y0){t= t0;y= y0
;}


/*:183*/
#line 181 "./solut.w"

}
#line 1 "./adabst.w"
/*:185*/

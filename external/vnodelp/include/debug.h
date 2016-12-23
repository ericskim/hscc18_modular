/*372:*/
#line 84 "./misc.w"

#ifndef DEBUG_H
#define DEBUG_H
#include "basiclinalg.h"
/*369:*/
#line 43 "./misc.w"

#ifdef VNODE_DEBUG
#define printMessage(s) {cerr << "\n *** " << __FILE__":"\
 <<__LINE__ << "   " << s << endl;}
#else
#define printMessage(s)  (0)
#endif

#ifdef VNODE_DEBUG
#define exitOnError(s) { printMessage(s); exit(-1);}
#else
#define exitOnError(s)   (0)
#endif


/*:369*/
#line 88 "./misc.w"


using namespace v_blas;
void checkIntersection(const iVector&a,const iVector&b);

#endif

/*:372*/

/*120:*/
#line 81 "./rounding.w"

#ifndef VNODEROUND_H
#define VNODEROUND_H

#ifdef FILIB_VNODE
#include <rounding_control/rounding_control_double.hpp>     
namespace v_bias{
/*117:*/
#line 32 "./rounding.w"


typedef filib::rounding_control<double,true> round_control;



/*:117*/
#line 88 "./rounding.w"

/*118:*/
#line 40 "./rounding.w"



inline void round_nearest()
{round_control::tonearest();}


inline void round_down()
{round_control::downward();}

inline void round_up()
{round_control::upward();}








/*:118*/
#line 89 "./rounding.w"
}
#endif


#ifdef PROFIL_VNODE
#include "Bias0.h"
namespace v_bias{
/*119:*/
#line 64 "./rounding.w"


inline void round_nearest()
{BiasRoundNear();}


inline void round_down()
{BiasRoundDown();}

inline void round_up()
{BiasRoundUp();}



/*:119*/
#line 96 "./rounding.w"

}
#endif



#endif
#line 97 "./vnode.w"



/*:120*/

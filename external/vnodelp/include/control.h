/*201:*/
#line 84 "./vnodecontrol.w"

#ifndef CONTROL_H
#define CONTROL_H
namespace vnodelp{
/*198:*/
#line 27 "./vnodecontrol.w"

typedef enum{first_entry,success,failure}Ind;


/*:198*/
#line 88 "./vnodecontrol.w"

/*199:*/
#line 43 "./vnodecontrol.w"
typedef enum{no,before_accept}
Interrupt;

/*:199*/
#line 89 "./vnodecontrol.w"

/*200:*/
#line 66 "./vnodecontrol.w"

class Control
{
public:

Ind ind;
Interrupt interrupt;
unsigned int order;
double atol,rtol;
double hmin;
Control():ind(first_entry),interrupt(no),
order(20),atol(1e-12),rtol(1e-12),
hmin(0){}
};

/*:200*/
#line 90 "./vnodecontrol.w"

}
#endif




#line 1 "./hoe.w"
/*:201*/

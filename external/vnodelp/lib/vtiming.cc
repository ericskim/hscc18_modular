/*376:*/
#line 17 "./timing.w"

#include <cassert> 
#include <sys/times.h> 
#include <unistd.h> 
#include <ctime> 

#include "vnodeinterval.h"
#include "vnoderound.h"

using namespace std;
static struct tms Tms;


double getTime()
{
times(&Tms);

v_bias::round_nearest();

long int ClockTcks= sysconf(_SC_CLK_TCK);
return(Tms.tms_utime)/(double(ClockTcks));
}

double getTotalTime(double start_time,double end_time)
{
assert(start_time<=end_time);
v_bias::round_nearest();

return(end_time-start_time);
}

#line 130 "./vnode.w"
/*:376*/

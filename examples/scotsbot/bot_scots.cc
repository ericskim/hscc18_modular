/*
 * mexReachController.cc
 *
 *  created on: 11.01.2016
 *      author: rungger
 */

//#include "mex.h"
#include "IO.hh"
//#include "ClassHandle.hh"
#include "UniformGrid.hh"
#include "Game.hh"
#include <khepera/khepera.h>

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netdb.h>
#include <khepera/khepera.h>
#include <signal.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netdb.h>
#include <sys/socket.h>


#define BUFLEN 100 // Buffer to receive UDP Data
#define SERVICE_PORT 5007 // Port for UDP Connection. Must be known to Client
#define UDP_PKT_SIZE 17


static knet_dev_t* dsPic; // robot pic microcontroller access
static int quitReq = 0; // quit variable for loop

using namespace std;

/*--------------------------------------------------------------------*/
/*!
 * Make sure the program terminate properly on a ctrl-c
 */
static void ctrlc_handler(int sig)
{
    quitReq = 1;

    kh4_set_speed(0, 0, dsPic); // stop robot
    kh4_SetMode(kh4RegIdle, dsPic);

    kh4_SetRGBLeds(0, 0, 0, 0, 0, 0, 0, 0, 0, dsPic); // clear rgb leds because consumes energy

    kb_change_term_mode(0); // revert to original terminal if called

    exit(0);
}


int main(void)

{


    typedef scots::Game Game;
    typedef scots::UniformGrid<std::vector<double> > Grid;

    const char* filename = "reach.scs";
    const char* filename2 = "reach2.scs";

    double kvr[3];
    
    int sw;
    
    sw=0;
    
    Game* game = new Game;
    Grid* ss = new Grid;;
    Grid* is = new Grid;;

    Game* gg = new Game;

    Game* gg2 = new Game;

    Grid* stateGrid = new Grid;
    Grid* inputGrid = new Grid;
    
    Grid* stateGrid2 = new Grid;
    Grid* inputGrid2 = new Grid;
    
    scots::IO::readControllerFromFile(gg, filename);
    scots::IO::readFromFile(stateGrid, filename, "state");
    scots::IO::readFromFile(inputGrid, filename, "input");

    scots::IO::readControllerFromFile(gg2, filename2);
    scots::IO::readFromFile(stateGrid2, filename2, "state");
    scots::IO::readFromFile(inputGrid2, filename2, "input");
     
    float vel_right, vel_left;

    int kp, ki, kd;
    int ctrloop_i=0;

    time_t start;


    int input = 0;
    char Buffer[100], revision, version;
    int int_t1 = 0;


    struct sockaddr_in myaddr, remaddr, myaddr1, remaddr1;
    int fd, fd1, i;
    unsigned int slen = sizeof(remaddr);
    char buf[BUFLEN]; /* message buffer */
    int recvlen; /* # bytes in acknowledgement message */
    char* server = "192.168.0.108"; /* change this to use a different server */
    double xp=0, y=0, theta=0;


    /* tuned parameters */
    kp = 10;
    ki = 5;
    kd = 1;
    kh4_ConfigurePID(kp, ki, kd, dsPic); // configure P,I,D

    printf("\nClosed Loop Trajectory Tracking Program\n");

    // initiate libkhepera and robot access
    if (kh4_init(0, NULL) != 0)
    {
        printf("\nERROR: could not initiate the libkhepera!\n\n");
        return -1;
    }

    /* open robot socket and store the handle in its pointer */
    dsPic = knet_open("Khepera4:dsPic", KNET_BUS_I2C, 0, NULL);

    if (dsPic == NULL)
    {
        printf("\nERROR: could not initiate communication with Kh4 dsPic\n\n");
        return -2;
    }

    // get revision
    if (kh4_revision(Buffer, dsPic) == 0)
    {
        version = (Buffer[0] >> 4) + 'A';
        revision = Buffer[0] & 0x0F;
        printf("\r\nVersion = %c, Revision = %u\r\n", version, revision);
    }

    signal(SIGINT, ctrlc_handler); // set signal for catching ctrl-c
    /*Pinging the server and getting x, y, theta values*/


    /* create a socket */

    if ((fd = socket(AF_INET, SOCK_DGRAM, 0)) == -1)
        printf("socket created\n");

    /* bind it to all local addresses and pick any port number */

    memset((char*)&myaddr, 0, sizeof(myaddr));
    myaddr.sin_family = AF_INET;
    myaddr.sin_addr.s_addr = htonl(INADDR_ANY);
    myaddr.sin_port = htons(SERVICE_PORT);

    struct timeval tv;
    tv.tv_sec = 3;
    tv.tv_usec = 0;
    if (setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv)) < 0)
    {
        perror("Error");
    }

    if (bind(fd, (struct sockaddr*)&myaddr, sizeof(myaddr)) < 0)
    {
        perror("bind failed");
        return 0;
    }

    /* now define remaddr, the address to whom we want to send messages */
    /* For convenience, the host address is expressed as a numeric IP address */
    /* that we will convert to a binary format via inet_aton */

    memset((char*)&remaddr, 0, sizeof(remaddr));
    remaddr.sin_family = AF_INET;
    remaddr.sin_port = htons(SERVICE_PORT);
    if (inet_aton(server, &remaddr.sin_addr) == 0)
    {
        fprintf(stderr, "inet_aton() failed\n");
        exit(1);
    }
    
   
for(int j=0;j<20;j++)

{
    if(sw==0)
    {

    game = gg;
    ss = stateGrid;
    is = inputGrid;
    } 

     if(sw==1)
    {

    game = gg2;
    ss = stateGrid2;
    is = inputGrid2;
    } 

    /* state space dimension */
    size_t sdim = ss->getDimension();
    size_t idim = is->getDimension();

    


    /* xyt should be found using UDP */


    /* now let's send the messages */
    for (;;)
    {
	printf("\nloop: %d\n", ctrloop_i++);

        recvlen = recvfrom(fd, buf, BUFLEN, 0, (struct sockaddr*)&remaddr, &slen);
        if (recvlen >= UDP_PKT_SIZE)
        {
            buf[recvlen] = 0;
            printf("\nReceived len: %d", recvlen);
            // printf("\nBuffer: %s", buf);
            sscanf(buf, "%lf %lf %lf", &xp, &y, &theta);
            // theta=theta*3.14/180;
            //xp = (xp + 121) / 10;
            //y = ((118 - y) / 10)-1.2;


		if(xp > 4 && xp < 6 && y > 4 && y < 6 && sw==0)
		{
		    printf("Target 1 Reached !!");
		    kh4_SetMode(kh4RegSpeed, dsPic);
		    kh4_set_speed(0, 0, dsPic);	
                    sw=1;
		    break;		
		}
             
               if(xp > 18.5 && xp < 20 && y > 16 && y < 18 && sw==1)
		{
		    printf("Target 2 Reached !!");
		    kh4_SetMode(kh4RegSpeed, dsPic);
		    kh4_set_speed(0, 0, dsPic);	
                    sw=0;
		    break;		
		}

        }
        else
        {

            printf("Stopping");
            kh4_SetMode(kh4RegSpeed, dsPic);
            kh4_set_speed(0, 0, dsPic);
        }


        /* current state*/

        kvr[0] = xp;
        kvr[1] = y;
        kvr[2] = theta;

        printf("\n%lf %lf %lf", xp, y, theta);

        std::vector<double> x(sdim);
        x.assign(kvr, kvr + 3);
        size_t xi;

        std::vector<size_t> ui;

        /* input */
        std::vector<double> u(idim);


        /* get index associated with  x */
        ss->xtoi(xi, x);

        /* get input idx associated with state index xi */
        ui = game->getLabel(xi);

        /* couldn't find input for this state*/
	if (recvlen >= UDP_PKT_SIZE)
	{
		if (ui.size() == 0)
		{
		    printf("\nERROR: could not find an input");

		    kh4_SetMode(kh4RegSpeed, dsPic);
		    kh4_set_speed(0, 0, dsPic);
		}
		else
		{ /* get input associated with index ui */
		    /* create matrix to store input */
		    printf("\nGO: got control agtion and sending to wheels");
		    for (size_t v = 0; v < ui.size(); v++)
		    {
		        is->itox(ui[v], u);
		    }

		    /* map u and v to vel_left and vel_right as follows */

		    vel_left = (((u[0] - 0.525 * u[1]) * 100) / KH4_SPEED_TO_MM_S);
		    vel_right = (((u[0] + 0.525 * u[1]) * 100) / KH4_SPEED_TO_MM_S);

		    printf("\n%f and %f\n", u[0], u[1]);
		}

	        kh4_SetMode(kh4RegSpeed, dsPic);
        	kh4_set_speed(vel_left, vel_right, dsPic);

	}
	else
	{
	    printf("\nno input will be sent as no paket reecieved !\n");		
	}



    }

printf("\nSwitching to the next target!\n");

}

    return 0;
}

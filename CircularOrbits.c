#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "philsplot.h"


/*
In this case:
u0 = R0/t0
u0 = sqrt(G*m0/x0)
x0 = R0
t0 = sqrt(xo^3/G*m0)
*/


//these are the working functions a that are being called
double starRHS (int i, double ux, double uy, double x, double y, double starx, double stary);                 //this evolves each star
double astRHS(int i, double ux, double uy, double x, double y, double arr[2][4]); //this evolves the asteroids
double initAsteroids(int n, double asteroid[4][500]);                             //this initializes the x,y,ux, and uy for the asteroids
                                                                                  //from here ux,uy,x,y will be referes to spatial components



int main(void){

  open_plot("1000x900")                ;//this is to open and draw axes in philsplot
  box_plot(-100., 100., -100., 100., 1., 6, "x", "y","","Evolution of binary star system with Asteroids");
  flush_plot();


  //Runge-Kutta vars
  double k_star[2][4][4];      // top row is for star 1, second row for star 2 and third for particle //in order, star #, spatials, RK component
  double k_ast[4][500][4];     // spatial down. asteroid number across, and K value in depth
  double step = .05;
  int K, s, o;                 // "s" is for spatial, "o" is for object, and "K" is the Kth value of k


  //position vars
  double star[2][4] = {{0. ,  0.5 ,  1. ,  0. },    //f[n][0] = ux, f[n][1] = uy, f[n][2] = x, f[n][3] = y
		      { 0. , -0.5 , -1. ,  0. }}; 


  double asteroid[4][500];    //for the nth asteroid, asteroid[0][n] = ux, asteroid[1][n] = uy, asteroid[2][n] = x, asteroid[3][n] = y. There are 500 asteroids, yea buddy
  int dead_ast[500];          //array for determining whether an asteroid was eaten/cast off or still orbiting
  int i;
  int numberAsteroids = 500;  //Again, there are 500 asteroids

  double t = 0;               //starting at time 0
  double x0;
  double u0;

  //calling initialized asteroids
  initAsteroids(numberAsteroids, asteroid);


  //Initializes dead asteroid array
  for (i = 0; i < numberAsteroids; i++){
    dead_ast[i] = 0;
  }



  //This is the meat of our program
    do {
     
      t = t + step;
 

     //BEGINNING OF STAR ROUTINE**********************************************************************************************


      for(K = 0; K < 4; K++){
	for(s = 0; s < 4; s++){
	  if (K == 0){
	    k_star[0][s][K] = step * starRHS(s, star[0][0], star[0][1], star[0][2], star[0][3], star[1][2], star[1][3]);
	    k_star[1][s][K] = step * starRHS(s, star[1][0], star[1][1], star[1][2], star[1][3], star[0][2], star[0][3]);
	  }
	  else if (K == 1){
	    k_star[0][s][K] = step * starRHS(s, star[0][0] + (0.5*k_star[0][s][K]), star[0][1] + (0.5*k_star[0][s][K]), star[0][2] + (0.5*k_star[0][s][K]), star[0][3] + (0.5*k_star[0][s][K]), star[1][2], star[1][3]);
	    k_star[1][s][K] = step * starRHS(s, star[1][0] + (0.5*k_star[1][s][K]), star[1][1] + (0.5*k_star[1][s][K]), star[1][2] + (0.5*k_star[1][s][K]), star[1][3] + (0.5*k_star[1][s][K]), star[0][2], star[0][3]);
	  }
	  else if (K == 2) {
	    k_star[0][s][K] = step * starRHS(s, star[0][0] + (0.5*k_star[0][s][K]), star[0][1] + (0.5*k_star[0][s][K]), star[0][2] + (0.5*k_star[0][s][K]), star[0][3] + (0.5*k_star[0][s][K]), star[1][2], star[1][3]);
	    k_star[1][s][K] = step * starRHS(s, star[1][0] + (0.5*k_star[1][s][K]), star[1][1] + (0.5*k_star[1][s][K]), star[1][2] + (0.5*k_star[1][s][K]), star[1][3] + (0.5*k_star[1][s][K]), star[0][2], star[0][3]);
	    }
	  else if (K == 3){
	    k_star[0][s][K] = step * starRHS(s, star[0][0] + k_star[0][s][K], star[0][1] + k_star[0][s][K], star[0][2] + k_star[0][s][K], star[0][3] + k_star[0][s][K], star[1][2], star[1][3]);
	    k_star[1][s][K] = step * starRHS(s, star[1][0] + k_star[1][s][K], star[1][1] + k_star[1][s][K], star[1][2] + k_star[1][s][K], star[1][3] + k_star[1][s][K], star[0][2], star[0][3]);

	    star[1][s] += (1./6.)*k_star[1][s][0] + (1./3.)*k_star[1][s][1] + (1./3.)*k_star[1][s][2] + (1./6.)*k_star[1][s][3];
	    star[0][s] += (1./6.)*k_star[0][s][0] + (1./3.)*k_star[0][s][1] + (1./3.)*k_star[0][s][2] + (1./6.)*k_star[0][s][3];
	  }
	}   //closes s loop
      }   //closes K loop



      putpoint_plot(star[0][2], star[0][3], 10, 1, 9., 0.2, 1); //draws star 1
      putpoint_plot(star[1][2], star[1][3], 10, 1, 9., 0.2, 1); //draws star 2     
      flush_plot(); //draws star 1 and 2 for real


     //END OF STAR ROUTINE***********************************************************************************************************************************



      //BEGINNING OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //evolving asteroids with RK routine
      for(K = 0; K < 4; K++){     //for Rk components
	for(s = 0; s < 4; s++){   //this is for the spatials
	  for(o = 0; o < numberAsteroids; o++){  //this is for each asteroid
	    if (dead_ast[o] == 1){ //checks to see if the asteroid is dead.by checking if it is either 1 or 0. if dead/ =1, then it moves on to the next one. 
	      continue;    //we learned a new command. this is it
	    }
	    else if (K == 0){
	      k_ast[s][o][K] = step * astRHS(s, asteroid[0][o], asteroid[1][o], asteroid[2][o], asteroid[3][o], star);
	    }
	    else if (K == 1){
	      k_ast[s][o][K] = step * astRHS(s, asteroid[0][o] + (0.5*k_ast[s][o][K]), asteroid[1][o] + (0.5*k_ast[s][o][K]), asteroid[2][o] + (0.5*k_ast[s][o][K]), asteroid[3][o] + (0.5*k_ast[s][o][K]), star);
	    }
	    else if (K == 2) {
	      k_ast[s][o][K] = step * astRHS(s, asteroid[0][o] + (0.5*k_ast[s][o][K]), asteroid[1][o] + (0.5*k_ast[s][o][K]), asteroid[2][o] + (0.5*k_ast[s][o][K]), asteroid[3][o] + (0.5*k_ast[s][o][K]), star);
	    }
	    else if (K == 3){
	      k_ast[s][o][K] = step *  astRHS(s, asteroid[0][o] + k_ast[s][o][K], asteroid[1][o] + k_ast[s][o][K], asteroid[2][o] + k_ast[s][o][K], asteroid[3][o] + k_ast[s][o][K], star);
	      asteroid[s][o] += (1./6.)*k_ast[s][o][0] + (1./3.)*k_ast[s][o][1] + (1./3.)*k_ast[s][o][2] + (1./6.)*k_ast[s][o][3];
	    }
	  }   //closes o loop
	}   //closes s loop
      }   //closes K loop


      for (o = 0; o < numberAsteroids;o++){
	putpoint_plot((2. * asteroid[2][o]),(2. * asteroid[3][o]), 1, 1, 2, 1, 1); //draws asteroid
	flush_plot();
      }

      //END OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////



      //This checks each asteroid to see if they are deadXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      double dist1, dist2; //this measures distance of ast from origin
      for (o = 0; o < (sizeof(asteroid[4])/sizeof(double)); o++){
	if(dead_ast[o] == 1){ //if the asteroid is already dead,move on
	  continue;
	}else{ //check to see if the asteroid is very far away from either star, if so, kill it.
	  dist1 = sqrt(pow((star[0][2]-asteroid[2][o]),2) + pow((star[0][3]-asteroid[3][o]),2));
	  dist2 = sqrt(pow((star[1][2]-asteroid[2][o]),2) + pow((star[1][3]-asteroid[3][o]),2));

	  if(dist1 > 750. && dist2 > 750.){ //if asteroid is too far away from system kill it
	    dead_ast[o] = 1;
	    printf("Asteroid %d has died from loneliness\n",o);
	  }
	  if(dist1 < 0.5 || dist2 < 0.5){ //if asteroid is too close to a star, the star plays hungry hungry hippo
	    dead_ast[o] = 1;
	    printf("Asteroid %d is dead. It was a tasty snack.\n",o);
	  }else{
	    continue;
	  }
	}
      }

      //Now Check to make sure asteroids are alive. If they're all dead, the program reinitializes the asteroid arry.//////////////
      int sum = 0;

      for (i = 0; i < numberAsteroids; i++){
	sum += dead_ast[i];
       }

      if (sum == numberAsteroids) {
	printf("All asteroids are dead, you monster\n");
	break;
      } else {
	sum = 0;
      }
      //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx

      //Now start erasing all the points that we drew. 
      //If this is not here the plot gets very crweded very quickly, albiet that it is spectacular getting there
      //doesnt erase every point for every t so trails that can be seen
      for(i=0;i<400;i++){
	delpoint_plot();
	flush_plot();
      }

    } while (t < 1000000000); //carry on forever

}//end main




//RHS for star
double starRHS(int i, double ux, double uy, double x, double y, double starx, double stary){
 if (i == 2){
  return ux;
 }
 if (i == 3){
   return uy;
 }
 if (i == 0){
   return (-(x - starx))/pow((pow(x-starx, 2) + pow(y-stary, 2)), 1.5);
 }
 if (i == 1){
   return (-(y - stary))/pow((pow(x-starx, 2) + pow(y-stary, 2)), 1.5);
 }
}



//RHS for asteroid
double astRHS(int i, double ux, double uy, double x, double y, double arr[2][4]) //arr[][] is referring to the star's components
{
 if (i == 2){
   return ux;
  }
  if (i == 3){
    return uy;
  }
  if (i == 0){
    return (-x)*((1/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 1.5)) + (1/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 1.5)));
  }
  if (i == 1){
    return (-y)*((1/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 1.5)) + (1/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 1.5)));
  }
} //end astRHS



//this initializes  every asteroids' spatials
double initAsteroids(int n, double asteroid[4][500]){
 int i, j;
 double d, theta;
 srand48(time(NULL)); //random number generators are for cool kids

 printf("Initializing asteroids....\n");


 for (i = 0; i < n; i++){            // cycles through each asteroid
   d = 250.*drand48();               //each asteroids is born within this radius from the origin
   theta = 2*M_PI*drand48(); 
   asteroid[2][i] = d*cos(theta);    // just setting random positions for x and y
   asteroid[3][i] = d*sin(theta);


   /*This segment now checks which quadrant the asteroid will land in and the velocites will be adjusted such that:
   Q1 --> ux negative, uy positive
   Q2 --> ux neg, uy neg
   Q3 --> ux pos, uy neg
   Q4 --> ux pos, uy pos
   and so that they are rotating in a counterclockwise fashion with the stars
   */

   if (asteroid[2][i] > 0 && asteroid[3][i] > 0){ //Q1
     asteroid[0][i] = -drand48() / 5.;
     asteroid[1][i] = drand48() / 5.;
   }
   else if(asteroid[2][i] < 0 && asteroid[3][i] > 0){ //Q2
     asteroid[0][i] = -drand48() / 5.;
     asteroid[1][i] = -drand48() / 5.;
   }
   else if(asteroid[2][i] < 0 && asteroid[3][i] < 0){ //Q3
     asteroid[0][i] = drand48() / 5.;
     asteroid[1][i] = -drand48() / 5.;
   }

   else if(asteroid[2][i] > 0 && asteroid[3][i] < 0){ //Q4
     asteroid[0][i] = drand48() / 5.;
     asteroid[1][i] = drand48()/ 5.;
   }
   else if(asteroid[2][i] == 0 || asteroid[3][i] == 0){ // for the lucky few to be born on an axis, they get ridiculous velocities. just because
     asteroid[0][i] = 9000*drand48();
     asteroid[1][i] = -9000*drand48();
   }

   printf("Asteroid %i spatial values are: %.3lf %.3lf %.3lf %.3lf\n", i, asteroid[0][i], asteroid[1][i], asteroid[2][i], asteroid[3][i] );
 } //end for loop
}//end InitAsteroids

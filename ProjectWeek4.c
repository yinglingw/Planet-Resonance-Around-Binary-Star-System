
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "philsplot.h"


/*
In this case: 
u0 = R0/t0 
u0 = sqrt(G*m0/x0)
x0 = R0
t0 = sqrt(xo^3/G*m0)
*/



double starRHS (int i, double ux, double uy, double x, double y);
double astRHS(int i, double ux, double uy, double x, double y, double arr[2][4]);
double initAsteroids(int n, double asteroid[4][200]);





int main(void){

  open_plot("1000x900");
  box_plot(-100., 100., -100., 100., 1., 6, "x", "y","","Evolution of binary star system with Asteroids");
  flush_plot();

  /*
  FILE *star1 = fopen("star1.txt", "w");
  FILE *star2 = fopen("star2.txt", "w");
  FILE *ASTEROID1 = fopen("asteroid1.txt", "w");
  FILE *ASTEROID2 = fopen("asteroid2.txt", "w");
  */
  /*
if (star1 == NULL || star2 == NULL || ASTEROID1 == NULL){
    printf("Can't open file\n");
    exit(1);
  }
  */
  
  // main vars
  double E = -0.5;
  double L = 1;
  


  //Runge-Kutta vars
  double k_star[2][4][4];   // Again,top row is for star 1, second row for star 2 and third for particle
  double k_ast[4][200][4]; // spatial down. asteroid number across, and K value in depth
  double step = .05;
  int K, s, o;    // "s" is for spatial, "o" is for object, and "K" is the Kth value of k


  //position vars
  double star[2][4] = {{0 ,  1.001 ,  1. ,  0 },   //f[n][0] = ux, f[n][1] = uy, f[n][2] = x, f[n][3] = y
		       { 0 , -1.001 , -1. ,  0 }};  //Top row is for star 1, second row for star 2 and third for particle;
  
  double asteroid[4][200]; //for the nth asteroid, asteroid[0][n] = ux, asteroid[1][n] = uy, asteroid[2][n] = x, asteroid[3][n] = y 
  int dead_ast[200];
  int i;
  int numberAsteroids = 200;

  double t = 0;
  double x0;
  double u0;
 
  
  initAsteroids(numberAsteroids, asteroid);




  //Initializes dead asteroid array/////////////////////////////////////////////////////////////////////////
  for (i = 0; i < numberAsteroids; i++){
    dead_ast[i] = 0;
  
    // printf("%d\n",dead_ast[i]);
	 }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////




     do {	
       
       t = t + step;
       


       //BEGINNING OF STAR ROUTINE**********************************************************************************************


       for(K = 0; K < 4; K++){
	      for(s = 0; s < 4; s++){
	              for(o = 0; o < 2; o++){
	            
			if (K == 0){
			  k_star[o][s][K] = step * starRHS(s, star[o][0], star[o][1], star[o][2], star[o][3]);
			} 
			else if (K == 1){
			  k_star[o][s][K] = step * starRHS(s, star[o][0] + (0.5*k_star[o][s][K]), star[o][1] + (0.5*k_star[o][s][K]), star[o][2] + (0.5*k_star[o][s][K]), star[o][3] + (0.5*k_star[o][s][K]));
			}
			else if (K == 2) {
			  k_star[o][s][K] = step * starRHS(s, star[o][0] + (0.5*k_star[o][s][K]), star[o][1] + (0.5*k_star[o][s][K]), star[o][2] + (0.5*k_star[o][s][K]), star[o][3] + (0.5*k_star[o][s][K]));
			} 
			else if (K == 3){
			  k_star[o][s][K] = step * starRHS(s, star[o][0] + k_star[o][s][K], star[o][1] + k_star[o][s][K], star[o][2] + k_star[o][s][K], star[o][3] + k_star[o][s][K]);
			  star[o][s] += (1./6.)*k_star[o][s][0] + (1./3.)*k_star[o][s][1] + (1./3.)*k_star[o][s][2] + (1./6.)*k_star[o][s][3];
			}
			
		      }   //closes o loop
	      }   //closes s loop
       }   //closes K loop

       
      
       putpoint_plot((2. * star[0][2]), (2. * star[0][3]), 1, 1, 5, 1.8, 1);
       delpoint_plot(); //draws star 1
       putpoint_plot((2. * star[1][2]), (2. * star[1][3]), 1, 1, 8, 0.5, 1); //draws star 2
       delpoint_plot();
       flush_plot();
       
       
       // fprintf(star1, "%lf %lf\n", star[0][2], star[0][3]);
       //  fprintf(star2, "%lf %lf\n", star[1][2], star[1][3]);

       //END OF STAR ROUTINE***********************************************************************************************************************************





       


       //BEGINNING OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       for(K = 0; K < 4; K++){
	      for(s = 0; s < 4; s++){
	              for(o = 0; o < numberAsteroids; o++){
		              if (dead_ast[o] == 1){
		         	continue;
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
	 putpoint_plot((2. * asteroid[2][o]),(2. * asteroid[3][o]), 1, 1, 2, 5., 5); //draws nth asteroid
       flush_plot();

       }
       //fprintf(ASTEROID1,"%lf %lf\n", asteroid[2][0], asteroid[3][0]);
       //fprintf(ASTEROID2,"%lf %lf\n", asteroid[2][1], asteroid[3][1]);
       //END OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////



       //This checks each asteroid to see if they are dead//////////////////////////////////////////////////////////////////////////
       double dist1, dist2; //this measures distance of ast from origin
       for (o = 0; o < (sizeof(asteroid[4])/sizeof(double)); o++){
	 if(dead_ast[o] == 1){
	   continue;
	 }else{
	 dist1 = sqrt(pow((star[0][2]-asteroid[2][o]),2) + pow((star[0][3]-asteroid[3][o]),2));
         dist2 = sqrt(pow((star[1][2]-asteroid[2][o]),2) + pow((star[1][3]-asteroid[3][o]),2));
 
	 if(dist1 > 200. && dist2 > 200.){ //if asteroid is too far away from system kill it
	   dead_ast[o] = 1;
	   printf("Asteroid %d has died from loneliness\n",o);
	 }
	 if(dist1 < 0.3 || dist2 < 0.3){ //if asteroid is too close to a star, the star plays hungry hungry hippo
	   dead_ast[o] = 1;
	   printf("Asteroid %d is dead. It was a tasty snack.\n",o);
	 }else{
	   continue;
	 }
	 }
       }

     //Checks to make sure all asteroids are alive. If they're all dead, the program reinitializes the asteroid arry.//////////////
       int sum = 0;

	for (i = 0; i < numberAsteroids; i++){
	  sum += dead_ast[i];
	  // printf("sum for i = %d is = %d\n",i,sum);
        }

	if (sum == numberAsteroids) {
	 printf("All asteroids are dead, you monster\n");
	 break;
       } else {
	 sum = 0;
       }
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      delay_plot(1);

     } while (t < 1000000000);


     //fclose(ASTEROID1);
     //fclose(ASTEROID2);
     //fclose(star1);
     //fclose(star2);


}

//RHS for star
double starRHS(int i, double ux, double uy, double x, double y){
   if (i == 2){
    return ux;
   }
   if (i == 3){
     return uy;
   }
   if (i == 0){
     return (-x)/pow((x*x + y*y), 1.5);
   }
   if (i == 1){
     return (-y)/pow((x*x + y*y), 1.5);
   }
}


//RHS for asteroid
double astRHS(int i, double ux, double uy, double x, double y, double arr[2][4])
{
  double c = 1.0; 
  if (i == 2){
    return ux;
   }
   if (i == 3){
     return uy;
   }
   if (i == 0){
     return (-c*x)*((1/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 3)) + (1/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 3)));
   }
   if (i == 1){
     return (-c*y)*((1/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 3)) + (1/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 3)));
   }
}


double initAsteroids(int n, double asteroid[4][200]){
  int i, j;
  double d, theta;
  srand48(time(NULL));

  printf("Initializing asteroids....\n");


  for (i = 0; i < n; i++){            // cycles through each asteroid
    d = 100.*drand48();
    theta = 2*M_PI*drand48();
    asteroid[2][i] = d*cos(theta);
    asteroid[3][i] = d*sin(theta);


    /*This segment checks which quadrant the asteroid will land in and the velocites will be adjusted such that:
    Q1 --> ux negative, uy positive
    Q2 --> ux neg, uy neg
    Q3 --> ux pos, uy neg
    Q4 --> ux pos, uy pos
    */


    if (asteroid[2][i] > 0 && asteroid[3][i] > 0){ //Q1
      asteroid[0][i] = -drand48() / 100.;
      asteroid[1][i] = drand48() / 100.;
    } 
    else if(asteroid[2][i] < 0 && asteroid[3][i] > 0){ //Q2
      asteroid[0][i] = -drand48() / 100.;
      asteroid[1][i] = -drand48() / 100.;
    }
    else if(asteroid[2][i] < 0 && asteroid[3][i] < 0){ //Q3
      asteroid[0][i] = drand48() / 100.;
      asteroid[1][i] = -drand48() / 100.;
    }

    else if(asteroid[2][i] > 0 && asteroid[3][i] < 0){ //Q4
      asteroid[0][i] = drand48() / 100.;
      asteroid[1][i] = drand48() / 100.;
    }
    else if(asteroid[2][i] == 0 || asteroid[3][i] == 0){
      asteroid[0][i] = 9000*drand48() / 100.;
      asteroid[1][i] = -9000*drand48() / 100.;
    }

    printf("Asteroid %i spatial values are: %.3lf %.3lf %.3lf %.3lf\n", i, asteroid[0][i], asteroid[1][i], asteroid[2][i], asteroid[3][i] );
  }
}

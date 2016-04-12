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



double starRHS (int i, double ux, double uy, double x, double y, double starx, double stary);
double astRHS(int i, double ux, double uy, double x, double y, double arr[2][4]);




int main(void){


FILE *converg1 = fopen("converge1.txt", "w");
FILE *converg2 = fopen("converge2.txt", "w");



if (converg1 == NULL || converg2 == NULL){
  printf("Can't open file\n");
  exit(1);
}



//Runge-Kutta vars
double k_star[2][4][4];   // Again,top row is for star 1, second row for star 2 and third for particle
double k_ast[4][2][4]; // spatial down. asteroid number across, and K value in depth
double step = .1;
int K, s, o;    // "s" is for spatial, "o" is for object, and "K" is the Kth value of k


//position vars
double star[2][4] = {{0 ,  .5 ,  1. ,  0 },   //f[n][0] = ux, f[n][1] = uy, f[n][2] = x, f[n][3] = y
		     {0 , -.5 , -1. ,  0 }};  //Top row is for star 1, second row for star 2 and third for particle;

double asteroid[4][500]; //for the nth asteroid, asteroid[0][n] = ux, asteroid[1][n] = uy, asteroid[2][n] = x, asteroid[3][n] = y
int dead_ast[500];
int i, j;
int numberAsteroids = 1;
int step_increase;


double E_star, E_ast;
double computedEast, computedEstar;
double errorAst, errorStar;
double t = 0;
double x0;
double u0;


 for(step_increase = 0; step_increase < 100001; step_increase++){
   
   
   asteroid[0][0] = -0.1;
   asteroid[1][0] = 0.;
   asteroid[2][0] = 0.;
   asteroid[3][0] = 10.;

   if (step_increase == 0){
       step = .01;
       printf("step intialized\n");
     } else {
     step = 6. * (pow(.071,(double)step_increase));
     }

   do {
   
     t += step;
   

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


     //END OF STAR ROUTINE***********************************************************************************************************************************
   

    

     //BEGINNING OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
     for(K = 0; K < 4; K++){
       for(s = 0; s < 4; s++){
	 for(o = 0; o < numberAsteroids; o++){
if (K == 0){
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

     //END OF ASTEROID ROUTINE////////////////////////////////////////////////////////////////////////////////////////////////////

     //printf("Can Runge-Kutta: %lf\n", asteroid[1][0] );
   } while(t < 1000.);


   if(step_increase == 0){ //fiducial value for star
     E_star = ((.5)*(sqrt((star[0][0]*star[0][0]) + (star[0][1]*star[0][1])))) + (1/sqrt((star[0][2]*star[0][2]) + (star[0][3]*star[0][3])));
     printf("Got fidcucial: %lf\n", E_star);
   } 
   else {
     computedEstar = ((.5)*(sqrt((star[0][0]*star[0][0]) + (star[0][1]*star[0][1])))) + (1/sqrt((star[0][2]*star[0][2]) + (star[0][3]*star[0][3])));
     errorStar = fabs(computedEstar - E_star)/E_star;
     
     // printf("Star:  %.16f %.16lf\n", errorStar, step);
     fprintf(converg1, "%.32lf %.32lf\n", step, errorStar);
   }



   if(step_increase == 0){ //fiducial value for asteroid
     E_ast = ((.5)*(sqrt((asteroid[0][0]*asteroid[0][0]) + (asteroid[1][0]*asteroid[1][0]))) + (1/sqrt((asteroid[2][0]*asteroid[2][0]) + (asteroid[3][0]*asteroid[3][0])) ));
     printf("Got fiducial: %lf\n", E_ast);
   } 
   else {
       computedEast = ((.5)*(sqrt((asteroid[0][0]*asteroid[0][0]) + (asteroid[1][0]*asteroid[1][0]))) + (1/sqrt((asteroid[2][0]*asteroid[2][0]) + (asteroid[3][0]*asteroid[3][0])) ));
     errorAst = fabs(computedEast - E_ast)/E_ast;
  
     // printf("Asteroid:  %.32lf %.32lf\n", errorAst, step);
     fprintf(converg2,"%.32lf %.32lf\n", step, errorAst);
     // printf ("Done with step # %d with h = %lf\n", step_increase, step);



   }

 
 }

     fclose(converg1);
     fclose(converg2);


}

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
   return (-(x - arr[1][2])/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 1.5)) + (-(x - arr[2][2])/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 1.5));
 }
 if (i == 1){
   return (-(y - arr[1][3])/pow(pow((x - arr[1][2]),2) + pow((y - arr[1][3]), 2), 1.5)) + (-(y - arr[2][3])/pow(pow((x - arr[2][2]),2) + pow((y - arr[2][3]), 2), 1.5));
 }
}


//periodic passive cross model, a[row][column]
//for 30X30 system we can go upto rho=0.14
//L(30),trelax(4),tmax(3),realz(3) = (10 data) 1 hr

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          50
# define trelax     10000
# define tmax       1000
# define realz      1000


int main()
{ clock_t tic = clock();

int    a[L][L],a_prime[L][L],b[L][L],position1[L*L]={0},position2[L*L]={0},NN,y,z,i,j,p,q,e;
long int t,tt,Tr,Tm,x1,x2;
double rho,c,d,slot1,slot2,slot3,slot4,block1,block2,rT=1.0,r,dt,pT,pA,pR,pTA,u,v,y1,y2,slope,vel,p0,pa,pb,pc;
char outputfile[200];
init_genrand(7687);
FILE   *fp;


rT=1.0;


fp = fopen("p0abc","w");
//sprintf(outputfile, "traj-%.2lf",rho);
//fp = fopen(outputfile,"w");


//loop for parameters starts
for(rho=0.01; rho<0.14; rho+=0.01){


NN = L*L*rho;
r = 4.0*NN*rT;
dt = 1.0/r;
Tr = (int)r*trelax;
Tm = (int)r*tmax;


p0=0; pa=0; pb=0; pc=0;


slot1 = 1.0/4.0;
slot2 = 2.0/4.0;
slot3 = 3.0/4.0;
slot4 = 4.0/4.0;



for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0;}
//random initialization of passive particles
  i=0; while (i<NN){
  y=genrand_real1()*L; z=genrand_real1()*L ;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){
  a[y][z]=2; a[(y+1)%L][z]=3; a[(y-1+L)%L][z]=3; a[y][(z+1)%L]=3; a[y][(z-1+L)%L]=3;
  position1[i] = y; position2[i] = z;
  i++ ; }
                  }


//realization loop starts
for(e=0;e<realz;e++){

//======================================================================================================================
//relaxation dynamics starts
for(t=0; t<Tr; t++){


//choosing rate to execute thermal or rotation move
d = genrand_real1();

//choosing particle for thermal move
i = genrand_real1()*NN;
p = position1[i]; q = position2[i];

  if(a[p][q]==2){

     //---------x+
     if(d < slot1){
       if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=2; 
                                         a[p][(q+2)%L]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q+1)%L]=3; 
                                         a[(p+1)%L][(q+1)%L]=3; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;      
            position1[i] = p; position2[i] = (q+1)%L;
                                         }  
                                      }

      //---------y+
      else if(d > slot1 && d < slot2){
        if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=2; 
                                         a[(p-2+L)%L][q]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q-1+L)%L]=3; 
                                         a[(p-1+L)%L][(q+1)%L]=3; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;       
                position1[i] = (p-1+L)%L; position2[i] = q;
                                         }  
                                      }
     //---------x-
     else if(d > slot2 && d < slot3){
       if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=2; 
                                         a[p][(q-2+L)%L]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q-1+L)%L]=3; 
                                         a[(p+1)%L][(q-1+L)%L]=3; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;             
            position1[i] = p; position2[i] = (q-1+L)%L;
                                         }  
                                      }
     //---------y-
     else if(d > slot3 && d < slot4){
       if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=2; 
                                         a[(p+2)%L][q]=3; 
                                         a[p][q]=3;
                                         a[(p+1)%L][(q-1+L)%L]=3; 
                                         a[(p+1)%L][(q+1)%L]=3; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                    
            position1[i] = (p+1)%L; position2[i] = q;   
                                         }  
                                      }

                          }//if choose thermal particle for thermal move

               }//Tr ends


//======================================================================================================================


u=0; v=0;

//measurement dynamics starts
for(t=0; t<Tm; t++){

//choosing rate to execute thermal or rotation move
d = genrand_real1();

//choosing particle for thermal move
i = genrand_real1()*NN;
p = position1[i]; q = position2[i];

  if(a[p][q]==2){

     //---------x+
     if(d < slot1){
       if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=2; 
                                         a[p][(q+2)%L]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q+1)%L]=3; 
                                         a[(p+1)%L][(q+1)%L]=3; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;      
            position1[i] = p; position2[i] = (q+1)%L;
                                         }  
                                      }

      //---------y+
      else if(d > slot1 && d < slot2){
        if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=2; 
                                         a[(p-2+L)%L][q]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q-1+L)%L]=3; 
                                         a[(p-1+L)%L][(q+1)%L]=3; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;       
                position1[i] = (p-1+L)%L; position2[i] = q;
                                         }  
                                      }
     //---------x-
     else if(d > slot2 && d < slot3){
       if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=2; 
                                         a[p][(q-2+L)%L]=3; 
                                         a[p][q]=3; 
                                         a[(p-1+L)%L][(q-1+L)%L]=3; 
                                         a[(p+1)%L][(q-1+L)%L]=3; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;             
            position1[i] = p; position2[i] = (q-1+L)%L;
                                         }  
                                      }
     //---------y-
     else if(d > slot3 && d < slot4){
       if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=2; 
                                         a[(p+2)%L][q]=3; 
                                         a[p][q]=3;
                                         a[(p+1)%L][(q-1+L)%L]=3; 
                                         a[(p+1)%L][(q+1)%L]=3; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                    
            position1[i] = (p+1)%L; position2[i] = q;   
                                         }  
                                      }

                          }//if choose thermal particle for thermal move


//measuring correlations after each mcs
i = 0; p = position1[i]; q = position2[i];

if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){p0++;}
if(a[(p-2+L)%L][(q+1)%L]==2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){pa++;}
if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]==2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){pb++;}
if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]==2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){pc++;}


            }//Tm ends

      }//realz ends 


//printing probabilities
fprintf(fp,"%lf  %lf  %lf  %lf  %lf\n",rho,p0/Tm/realz,pa/Tm/realz,pb/Tm/realz,pc/Tm/realz);

         }//parameter loop ends


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends

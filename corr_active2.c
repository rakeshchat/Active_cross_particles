//periodic active-cross model with 1 active particle and activity in the positive X-direction, a[row][column]
//for 30X30 system we can go upto rho=0.14
//L(30),trelax(4),tmax(3),realz(3) = (10 data) 1 hr

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          100
# define trelax     10000
# define tmax       10000
# define realz      20


int main()
{ clock_t tic = clock();

int    a[L][L],a_prime[L][L],b[L][L],position1[L*L]={0},position2[L*L]={0},NN,y,z,i,j,p,q,e;
long int t,tt,Tr,Tm,x1,x2;
double rho,c,d,slot1,slot2,slot3,slot4,rT=1.0,rA,rR,r,dt,pT,pA,pR,pTA,u,v,y1,y2,slope,vel,C_plus,C_minus;
double Cp_a,Cp_b,Cp_c,Cp_d,Cp_e,Cp_ac,Cp_ad,Cp_ae,Cp_be,Cp_ce,Cp_ace,Cm_a,Cm_b,Cm_c,Cm_d,Cm_e,Cm_ac,Cm_ad,Cm_ae,Cm_be,Cm_ce,Cm_ace;
double Cp_a1,Cp_b1,Cp_c1,Cp_d1,Cp_e1,Cp_ac1,Cp_ad1,Cp_ae1,Cp_be1,Cp_ce1,Cp_ace1,Cm_a1,Cm_b1,Cm_c1,Cm_d1,Cm_e1,Cm_ac1,Cm_ad1,Cm_ae1,Cm_be1,Cm_ce1,Cm_ace1;
char outputfile[200];
init_genrand(76879);
FILE   *fp;


rho=0.130;
//rA=10.0;
rR=0.0;


sprintf(outputfile, "corrdata-rA-%.3lfa",rho);
fp = fopen(outputfile,"w");


//loop for parameters starts
for(rA=0.01; rA<10.0; rA*=1.5){
//for(rho=0.01; rho<=0.13; rho+=0.01){


NN = L*L*rho;
r = 4.0*NN*rT + rA + 2.0*rR;
dt = 1.0/r;
Tr = (int)r*trelax;
Tm = (int)r*tmax;


vel=0; C_plus=0; C_minus=0;


Cp_a=0; Cp_b=0; Cp_c=0; Cp_d=0; Cp_e=0;
Cp_ac=0; Cp_ad=0; Cp_ae=0; Cp_be=0; Cp_ce=0; Cp_ace=0; 
Cm_a=0; Cm_b=0; Cm_c=0; Cm_d=0; Cm_e=0;
Cm_ac=0; Cm_ad=0; Cm_ae=0; Cm_be=0; Cm_ce=0; Cm_ace=0;
Cp_a1=0; Cp_b1=0; Cp_c1=0; Cp_d1=0; Cp_e1=0;
Cp_ac1=0; Cp_ad1=0; Cp_ae1=0; Cp_be1=0; Cp_ce1=0; Cp_ace1=0; 
Cm_a1=0; Cm_b1=0; Cm_c1=0; Cm_d1=0; Cm_e1=0;
Cm_ac1=0; Cm_ad1=0; Cm_ae1=0; Cm_be1=0; Cm_ce1=0; Cm_ace1=0;


pT = 4.0*NN*rT/r;
pA = 1.0*rA/r;
pR = 2.0*rR/r;

pTA = pT + pA;

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

//randomly choose one active among all passive particles
  i=0; while (i<1){
  y = position1[i]; z = position2[i];
  a[y][z]=5; a[(y+1)%L][z]=6; a[(y-1+L)%L][z]=6; a[y][(z+1)%L]=6; a[y][(z-1+L)%L]=6;
  i++ ;   
                  }

//realization loop starts
for(e=0;e<realz;e++){

//======================================================================================================================
//relaxation dynamics starts
for(t=0; t<Tr; t++){


//choosing which move will take place among thermal,active or rotation move
c = genrand_real1();
//choosing rate to execute thermal or rotation move
d = genrand_real1();

//for thermal move
if(c < pT){
//choosing particle for thermal move
i = genrand_real1()*NN;
p = position1[i]; q = position2[i];

  if(a[p][q]==2){

     //---------x+
     if(d < slot1){
       if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 && 
          a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 && 
          a[(p-2+L)%L][(q+1)%L]!=7 && a[(p-1+L)%L][(q+2)%L]!=7 && a[p][(q+3)%L]!=7 && a[(p+1)%L][(q+2)%L]!=7 && a[(p+2)%L][(q+1)%L]!=7 && 
          a[(p-2+L)%L][(q+1)%L]!=8 && a[(p-1+L)%L][(q+2)%L]!=8 && a[p][(q+3)%L]!=8 && a[(p+1)%L][(q+2)%L]!=8 && a[(p+2)%L][(q+1)%L]!=8 && 
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9){
             
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
        if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 && 
           a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 && 
           a[(p-2+L)%L][(q+1)%L]!=7 && a[(p-1+L)%L][(q+2)%L]!=7 && a[(p-3+L)%L][q]!=7 && a[(p-2+L)%L][(q-1+L)%L]!=7 && a[(p-1+L)%L][(q-2+L)%L]!=7 && 
           a[(p-2+L)%L][(q+1)%L]!=8 && a[(p-1+L)%L][(q+2)%L]!=8 && a[(p-3+L)%L][q]!=8 && a[(p-2+L)%L][(q-1+L)%L]!=8 && a[(p-1+L)%L][(q-2+L)%L]!=8 && 
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9){
             
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
       if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && 
          a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && 
          a[(p-1+L)%L][(q-2+L)%L]!=7 && a[(p-2+L)%L][(q-1+L)%L]!=7 && a[p][(q-3+L)%L]!=7 && a[(p+1)%L][(q-2+L)%L]!=7 && a[(p+2)%L][(q-1+L)%L]!=7 && 
          a[(p-1+L)%L][(q-2+L)%L]!=8 && a[(p-2+L)%L][(q-1+L)%L]!=8 && a[p][(q-3+L)%L]!=8 && a[(p+1)%L][(q-2+L)%L]!=8 && a[(p+2)%L][(q-1+L)%L]!=8 && 
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9){
             
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
       if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && 
          a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && 
          a[(p+2)%L][(q+1)%L]!=7 && a[(p+1)%L][(q+2)%L]!=7 && a[(p+3)%L][q]!=7 && a[(p+2)%L][(q-1+L)%L]!=7 && a[(p+1)%L][(q-2+L)%L]!=7 && 
          a[(p+2)%L][(q+1)%L]!=8 && a[(p+1)%L][(q+2)%L]!=8 && a[(p+3)%L][q]!=8 && a[(p+2)%L][(q-1+L)%L]!=8 && a[(p+1)%L][(q-2+L)%L]!=8 && 
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9){
             
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

     else if(a[p][q]==5){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=5; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;
            position1[i] = p; position2[i] = (q+1)%L;                        
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=5; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                 
            position1[i] = (p-1+L)%L; position2[i] = q; 
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=5; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
            position1[i] = p; position2[i] = (q-1+L)%L;          
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=5; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                        
            position1[i] = (p+1)%L; position2[i] = q;
                                         }  
                                      }

                           }//if choose active x+ particle for thermal move ends

      else if(a[p][q]==7){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=7; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                    
            position1[i] = p; position2[i] = (q+1)%L;           
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=7; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;    
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=7; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
                position1[i] = p; position2[i] = (q-1+L)%L;     
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=7; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
                position1[i] = (p+1)%L; position2[i] = q;      
                                         }  
                                      }

                           }//if choose active y+ particle for thermal move ends

      else if(a[p][q]==8){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=8; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                    
                position1[i] = p; position2[i] = (q+1)%L;      
                          
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=8; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;     
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=8; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;
                position1[i] = p; position2[i] = (q-1+L)%L;                                       
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=8; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                    
                position1[i] = (p+1)%L; position2[i] = q;                             
                                         }  
                                      }

                           }//if choose active x- particle for thermal move ends


      else if(a[p][q]==9){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=9; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                      
                position1[i] = p; position2[i] = (q+1)%L;     
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=9; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;   
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=9; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
                position1[i] = p; position2[i] = (q-1+L)%L;                     
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=9; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
                position1[i] = (p+1)%L; position2[i] = q;                               
                                         }  
                                      }

                           }//if choose active y- particle for thermal move ends

                     }//total thermal move ends


//for active moves
else if(c > pT && c < pTA){
//choosing active particle for active move
i = 0;
p = position1[i]; q = position2[i];

   //active move towards x+ 
   if(a[p][q]==5){
            if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
                                         a[p][(q+1)%L]=5; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                      
            position1[i] = p; position2[i] = (q+1)%L;                                                      
                                         }       
                          
                 }

    //active move towards y+ 
    else if(a[p][q]==7){
            if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
                                         a[(p-1+L)%L][q]=7; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
            position1[i] = (p-1+L)%L; position2[i] = q;                         
                                         }       
                          
                     }

    //active move towards x-
    else if(a[p][q]==8){
            if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
                                         a[p][(q-1+L)%L]=8; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
            position1[i] = p; position2[i] = (q-1+L)%L;                               
                                         }     
                          
                     }

    //active move towards y-
    else if(a[p][q]==9){
            if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
                                         a[(p+1)%L][q]=9; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
            position1[i] = (p+1)%L; position2[i] = q;                                      
                                         }      
                          
                     }

             }//total active move ends

//for rotation moves
else if(c > pTA){
//choosing active particle for active move
i = 0;
p = position1[i]; q = position2[i];

      //rotational move of active particle with x+
      if(a[p][q]==5){
        if(d < 0.5){a[p][q]=7;  }
        else {a[p][q]=9;  }
              }

      //rotational move of active particle with y+
      else if(a[p][q]==7){
        if(d < 0.5){a[p][q]=5; }
        else {a[p][q]=8; }
              }

     //rotational move of active particle with x-
     else if(a[p][q]==8){
        if(d < 0.5){a[p][q]=7; }
        else {a[p][q]=9; }
              }

     //rotational move of active particle with y-
     else if(a[p][q]==9){
        if(d < 0.5){a[p][q]=5; }
        else {a[p][q]=8; }
              }

         }//total rotational move ends


     }//Tr ends


//for(y=0;y<L;y++)for(z=0;z<L;z++){a_prime[y][z] = a[y][z];}
//for(y=0;y<L;y++){for(z=0;z<L;z++){printf("%d",a_prime[y][z]); } printf("\n"); }

//======================================================================================================================

//for(y=0;y<L;y++){for(z=0;z<L;z++){printf("%d",a[y][z]); } printf("\n"); } printf("\n"); 
//for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z] = 0;}
//for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z] = a_prime[y][z];}


u=0; v=0;

//measurement dynamics starts
for(t=0; t<Tm; t++){

//choosing which move will take place among thermal,active or rotation move
c = genrand_real1();
//choosing rate to execute thermal or rotation move
d = genrand_real1();

//for thermal move
if(c < pT){
//choosing particle for thermal move
i = genrand_real1()*NN;
p = position1[i]; q = position2[i];

  if(a[p][q]==2){

     //---------x+
     if(d < slot1){
       if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 && 
          a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 && 
          a[(p-2+L)%L][(q+1)%L]!=7 && a[(p-1+L)%L][(q+2)%L]!=7 && a[p][(q+3)%L]!=7 && a[(p+1)%L][(q+2)%L]!=7 && a[(p+2)%L][(q+1)%L]!=7 && 
          a[(p-2+L)%L][(q+1)%L]!=8 && a[(p-1+L)%L][(q+2)%L]!=8 && a[p][(q+3)%L]!=8 && a[(p+1)%L][(q+2)%L]!=8 && a[(p+2)%L][(q+1)%L]!=8 && 
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9){
             
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
        if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 && 
           a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 && 
           a[(p-2+L)%L][(q+1)%L]!=7 && a[(p-1+L)%L][(q+2)%L]!=7 && a[(p-3+L)%L][q]!=7 && a[(p-2+L)%L][(q-1+L)%L]!=7 && a[(p-1+L)%L][(q-2+L)%L]!=7 && 
           a[(p-2+L)%L][(q+1)%L]!=8 && a[(p-1+L)%L][(q+2)%L]!=8 && a[(p-3+L)%L][q]!=8 && a[(p-2+L)%L][(q-1+L)%L]!=8 && a[(p-1+L)%L][(q-2+L)%L]!=8 && 
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9){
             
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
       if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && 
          a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && 
          a[(p-1+L)%L][(q-2+L)%L]!=7 && a[(p-2+L)%L][(q-1+L)%L]!=7 && a[p][(q-3+L)%L]!=7 && a[(p+1)%L][(q-2+L)%L]!=7 && a[(p+2)%L][(q-1+L)%L]!=7 && 
          a[(p-1+L)%L][(q-2+L)%L]!=8 && a[(p-2+L)%L][(q-1+L)%L]!=8 && a[p][(q-3+L)%L]!=8 && a[(p+1)%L][(q-2+L)%L]!=8 && a[(p+2)%L][(q-1+L)%L]!=8 && 
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9){
             
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
       if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && 
          a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && 
          a[(p+2)%L][(q+1)%L]!=7 && a[(p+1)%L][(q+2)%L]!=7 && a[(p+3)%L][q]!=7 && a[(p+2)%L][(q-1+L)%L]!=7 && a[(p+1)%L][(q-2+L)%L]!=7 && 
          a[(p+2)%L][(q+1)%L]!=8 && a[(p+1)%L][(q+2)%L]!=8 && a[(p+3)%L][q]!=8 && a[(p+2)%L][(q-1+L)%L]!=8 && a[(p+1)%L][(q-2+L)%L]!=8 && 
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9){
             
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

     else if(a[p][q]==5){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=5; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;
            position1[i] = p; position2[i] = (q+1)%L;          u++;                
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=5; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;    v++;             
            position1[i] = (p-1+L)%L; position2[i] = q; 
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=5; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
            position1[i] = p; position2[i] = (q-1+L)%L;          u--;
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=5; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;         v--;               
            position1[i] = (p+1)%L; position2[i] = q;
                                         }  
                                      }

                           }//if choose active x+ particle for thermal move ends

      else if(a[p][q]==7){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=7; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                    
            position1[i] = p; position2[i] = (q+1)%L;           u++;
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=7; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;    v++;
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=7; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
                position1[i] = p; position2[i] = (q-1+L)%L;     u--;
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=7; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
                position1[i] = (p+1)%L; position2[i] = q;      v--;
                                         }  
                                      }

                           }//if choose active y+ particle for thermal move ends

      else if(a[p][q]==8){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=8; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                    
                position1[i] = p; position2[i] = (q+1)%L;      u++;
                          
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=8; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;     v++;
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=8; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;
                position1[i] = p; position2[i] = (q-1+L)%L;    u--;                                    
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=8; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                    
                position1[i] = (p+1)%L; position2[i] = q;    v--;                         
                                         }  
                                      }

                           }//if choose active x- particle for thermal move ends


      else if(a[p][q]==9){

        //---------x+
        if(d < slot1){
          if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
             
                                         a[p][(q+1)%L]=9; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                      
                position1[i] = p; position2[i] = (q+1)%L;     u++;
                                         }  
                                      }
        //---------y+
        else if(d > slot1 && d < slot2){
               if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
             
                                         a[(p-1+L)%L][q]=9; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
                position1[i] = (p-1+L)%L; position2[i] = q;   v++;
                                         }  
                                      }
        //---------x-
        else if(d > slot2 && d < slot3){
               if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
             
                                         a[p][(q-1+L)%L]=9; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
                position1[i] = p; position2[i] = (q-1+L)%L;     u--;                 
                                         }  
                                      }
        //---------y-
        else if(d > slot3 && d < slot4){
               if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
             
                                         a[(p+1)%L][q]=9; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
                position1[i] = (p+1)%L; position2[i] = q;     v--;                          
                                         }  
                                      }

                           }//if choose active y- particle for thermal move ends

                     }//total thermal move ends


//for active moves
else if(c > pT && c < pTA){
//choosing active particle for active move
i = 0;
p = position1[i]; q = position2[i];

   //active move towards x+ 
   if(a[p][q]==5){
            if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2){
                                         a[p][(q+1)%L]=5; 
                                         a[p][(q+2)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                      
            position1[i] = p; position2[i] = (q+1)%L;          u++;                                            
                                         }       
                          
                 }

    //active move towards y+ 
    else if(a[p][q]==7){
            if(a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2){
                                         a[(p-1+L)%L][q]=7; 
                                         a[(p-2+L)%L][q]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][(q+1)%L]=6; 
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;                     
            position1[i] = (p-1+L)%L; position2[i] = q;        v++;                 
                                         }       
                          
                     }

    //active move towards x-
    else if(a[p][q]==8){
            if(a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2){
                                         a[p][(q-1+L)%L]=8; 
                                         a[p][(q-2+L)%L]=6; 
                                         a[p][q]=6; 
                                         a[(p-1+L)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;                      
            position1[i] = p; position2[i] = (q-1+L)%L;         u--;                          
                                         }     
                          
                     }

    //active move towards y-
    else if(a[p][q]==9){
            if(a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2){
                                         a[(p+1)%L][q]=9; 
                                         a[(p+2)%L][q]=6; 
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6; 
                                         a[(p+1)%L][(q+1)%L]=6; 
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;                       
            position1[i] = (p+1)%L; position2[i] = q;          v--;                             
                                         }      
                          
                     }

             }//total active move ends

//for rotation moves
else if(c > pTA){
//choosing active particle for active move
i = 0;
p = position1[i]; q = position2[i];

      //rotational move of active particle with x+
      if(a[p][q]==5){
        if(d < 0.5){a[p][q]=7;  }
        else {a[p][q]=9;  }
              }

      //rotational move of active particle with y+
      else if(a[p][q]==7){
        if(d < 0.5){a[p][q]=5; }
        else {a[p][q]=8; }
              }

     //rotational move of active particle with x-
     else if(a[p][q]==8){
        if(d < 0.5){a[p][q]=7; }
        else {a[p][q]=9; }
                   }

     //rotational move of active particle with y-
     else if(a[p][q]==9){
        if(d < 0.5){a[p][q]=5; }
        else {a[p][q]=8; }
                   }

                }//total rotational move ends


//measuring correlations after each mcs
for(y=0;y<L;y++){for(z=0;z<L;z++){if(a[y][z]==5){p=y; q=z;} } }

for(y=0;y<L;y++){for(z=0;z<L;z++){b[y][z] = a[y][z];} }


if(b[(p-2+L)%L][(q+1)%L]==2){b[(p-2+L)%L][(q+1)%L]=1;}else{b[(p-2+L)%L][(q+1)%L]=0;}
if(b[(p-1+L)%L][(q+2)%L]==2){b[(p-1+L)%L][(q+2)%L]=1;}else{b[(p-1+L)%L][(q+2)%L]=0;}
if(b[p][(q+3)%L]==2){b[p][(q+3)%L]=1;}else{b[p][(q+3)%L]=0;}
if(b[(p+1)%L][(q+2)%L]==2){b[(p+1)%L][(q+2)%L]=1;}else{b[(p+1)%L][(q+2)%L]=0;}
if(b[(p+2)%L][(q+1)%L]==2){b[(p+2)%L][(q+1)%L]=1;}else{b[(p+2)%L][(q+1)%L]=0;}


if(b[(p-2+L)%L][(q-1+L)%L]==2){b[(p-2+L)%L][(q-1+L)%L]=1;}else{b[(p-2+L)%L][(q-1+L)%L]=0;}
if(b[(p-1+L)%L][(q-2+L)%L]==2){b[(p-1+L)%L][(q-2+L)%L]=1;}else{b[(p-1+L)%L][(q-2+L)%L]=0;}
if(b[p][(q-3+L)%L]==2){b[p][(q-3+L)%L]=1;}else{b[p][(q-3+L)%L]=0;}
if(b[(p+1)%L][(q-2+L)%L]==2){b[(p+1)%L][(q-2+L)%L]=1;}else{b[(p+1)%L][(q-2+L)%L]=0;}
if(b[(p+2)%L][(q-1+L)%L]==2){b[(p+2)%L][(q-1+L)%L]=1;}else{b[(p+2)%L][(q-1+L)%L]=0;}


Cp_a += b[(p-2+L)%L][(q+1)%L]; 
Cp_b += b[(p-1+L)%L][(q+2)%L]; 
Cp_c += b[p][(q+3)%L]; 
Cp_d += b[(p+1)%L][(q+2)%L]; 
Cp_e += b[(p+2)%L][(q+1)%L];

Cp_ac += b[(p-2+L)%L][(q+1)%L]*b[p][(q+3)%L];
Cp_ad += b[(p-2+L)%L][(q+1)%L]*b[(p+1)%L][(q+2)%L];
Cp_ae += b[(p-2+L)%L][(q+1)%L]*b[(p+2)%L][(q+1)%L];
Cp_be += b[(p-1+L)%L][(q+2)%L]*b[(p+2)%L][(q+1)%L];
Cp_ce += b[p][(q+3)%L]*b[(p+2)%L][(q+1)%L];
Cp_ace += b[(p-2+L)%L][(q+1)%L]*b[p][(q+3)%L]*b[(p+2)%L][(q+1)%L];


Cm_a += b[(p-2+L)%L][(q-1+L)%L]; 
Cm_b += b[(p-1+L)%L][(q-2+L)%L]; 
Cm_c += b[p][(q-3+L)%L]; 
Cm_d += b[(p+1)%L][(q-2+L)%L]; 
Cm_e += b[(p+2)%L][(q-1+L)%L];

Cm_ac += b[(p-2+L)%L][(q-1+L)%L]*b[p][(q-3+L)%L];
Cm_ad += b[(p-2+L)%L][(q-1+L)%L]*b[(p+1)%L][(q-2+L)%L];
Cm_ae += b[(p-2+L)%L][(q-1+L)%L]*b[(p+2)%L][(q-1+L)%L];
Cm_be += b[(p-1+L)%L][(q-2+L)%L]*b[(p+2)%L][(q-1+L)%L];
Cm_ce += b[p][(q-3+L)%L]*b[(p+2)%L][(q-1+L)%L];
Cm_ace += b[(p-2+L)%L][(q-1+L)%L]*b[p][(q-3+L)%L]*b[(p+2)%L][(q-1+L)%L];


            }//Tm ends

      }//realz ends


Cp_a1 = Cp_a/1.0/Tm/realz; 
Cp_b1 = Cp_b/1.0/Tm/realz; 
Cp_c1 = Cp_c/1.0/Tm/realz; 
Cp_d1 = Cp_d/1.0/Tm/realz; 
Cp_e1 = Cp_e/1.0/Tm/realz;
Cp_ac1 = Cp_ac/1.0/Tm/realz; 
Cp_ad1 = Cp_ad/1.0/Tm/realz; 
Cp_ae1 = Cp_ae/1.0/Tm/realz; 
Cp_be1 = Cp_be/1.0/Tm/realz; 
Cp_ce1 = Cp_ce/1.0/Tm/realz; 
Cp_ace1 = Cp_ace/1.0/Tm/realz;


Cm_a1 = Cm_a/1.0/Tm/realz; 
Cm_b1 = Cm_b/1.0/Tm/realz; 
Cm_c1 = Cm_c/1.0/Tm/realz; 
Cm_d1 = Cm_d/1.0/Tm/realz; 
Cm_e1 = Cm_e/1.0/Tm/realz;
Cm_ac1 = Cm_ac/1.0/Tm/realz; 
Cm_ad1 = Cm_ad/1.0/Tm/realz; 
Cm_ae1 = Cm_ae/1.0/Tm/realz; 
Cm_be1 = Cm_be/1.0/Tm/realz; 
Cm_ce1 = Cm_ce/1.0/Tm/realz; 
Cm_ace1 = Cm_ace/1.0/Tm/realz;




C_plus = (1.0 - Cp_a1 - Cp_b1 - Cp_c1 - Cp_d1 - Cp_e1 + Cp_ac1 + Cp_ad1 + Cp_ae1 + Cp_be1 + Cp_ce1 - Cp_ace1);
C_minus = (1.0 - Cm_a1 - Cm_b1 - Cm_c1 - Cm_d1 - Cm_e1 + Cm_ac1 + Cm_ad1 + Cm_ae1 + Cm_be1 + Cm_ce1 - Cm_ace1);
vel = 1.0*(rA + rT)*C_plus - 1.0*rT*C_minus;



//printing correlation
//fprintf(fp,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",Cp_a1,Cp_b1,Cp_c1,Cp_d1,Cp_e1,Cp_ac1,Cp_ad1,Cp_ae1,Cp_be1,Cp_ce1,Cp_ace1,Cm_a1,Cm_b1,Cm_c1,Cm_d1,Cm_e1,Cm_ac1,Cm_ad1,Cm_ae1,Cm_be1,Cm_ce1,Cm_ace1,rA);


//printing velocity from correlation
fprintf(fp,"%lf   %lf    %lf    %lf\n",rA,C_plus,C_minus,vel);

    }//parameter loop ends


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends

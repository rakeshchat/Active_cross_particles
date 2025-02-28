//periodic active-cross model with 1 active particle and activity in the positive X-direction, a[row][column]
//for 30X30 system we can go upto rho=0.14
//L(30),trelax(4),tmax(3),realz(2) = gcc(4 min) icc(1 min 40 sec)

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

int    a[L][L],a_prime[L][L],position1[L*L]={0},position2[L*L]={0},NN,y,z,i,j,p,q,e;
long int t,tt,Tr,Tm,x1,x2;
double rho,c,d,slot1,slot2,slot3,slot4,block1,block2,rT,rA,rR,r,dt,pT,pA,pR,pTA,u,v,y1,y2,slope;
double *s,*xt,*yt,*xtav,*ytav,*msd,*msdav;
char outputfile[200];
init_genrand(7687);
FILE   *fp;


rho=0.06;
rT=1.0;


fp = fopen("asd-psv","w");
//sprintf(outputfile, "traj-%.2lf",rho);
//fp = fopen(outputfile,"w");


//loop for parameters starts
//for(rA=0.20; rA<30.0; rA*=1.3){
//for(rho=0.01; rho<0.12; rho+=0.05){



NN = L*L*rho;
r = 4.0*NN*rT;
dt = 1.0/r;
Tr = (int)r*trelax;
Tm = (int)r*tmax;


s=(double *)calloc(Tm,sizeof(double));
xt=(double *)calloc(Tm,sizeof(double));
yt=(double *)calloc(Tm,sizeof(double));
xtav=(double *)calloc(Tm,sizeof(double));
ytav=(double *)calloc(Tm,sizeof(double));
msd=(double *)calloc(Tm,sizeof(double));
msdav=(double *)calloc(Tm,sizeof(double));
for(t=0;t<Tm;t++){s[t]=0; xt[t]=0; yt[t]=0; xtav[t]=0; ytav[t]=0; msd[t]=0; msdav[t]=0;}
slope=0;


pT = 4.0*NN*rT/r;


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

//randomly choose one among all passive particles
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

       xt[t] = u; yt[t] = v;

       xtav[t] += xt[t]*xt[t];  ytav[t] += yt[t]*yt[t];

            }//Tm ends

      }//realz ends


//printing MSD for rotating active tracer
for(t=0;i<Tm;t++){fprintf(fp,"%lf    %lf\n",1.0*t*dt,(xtav[t] + ytav[t])/realz);}


//         }//parameter loop ends


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
free(s); free(xt); free(yt); free(xtav); free(ytav); free(msd); free(msdav);
}//main ends

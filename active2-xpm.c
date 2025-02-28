//                 _____
//                |     |
//                |  3  |
//            ____|     |____
//           |               |
//           | 3     2     3 |        passive
//           |____       ____|
//                |     |
//                |  3  |
//                |_____|      

//                 _____
//                |     |
//                |  6  |
//            ____|     |____
//           |               |
//           | 6     5     6 | -->>   active in x+
//           |____       ____|
//                |     |
//                |  6  |
//                |_____|    


//                   ^
//                   |
//                 __|__
//                |     |
//                |  6  |
//            ____|     |____
//           |               |
//           | 6     7     6 |        active in y+
//           |____       ____|
//                |     |
//                |  6  |
//                |_____|    


//                 _____
//                |     |
//                |  6  |
//            ____|     |____
//           |               |
//      <<-- | 6     8     6 |        active in x-
//           |____       ____|
//                |     |
//                |  6  |
//                |_____|   


//                 _____
//                |     |
//                |  6  |
//            ____|     |____
//           |               |
//           | 6     9     6 |        active in y-
//           |____       ____|
//                |     |
//                |  6  |
//                |_____|   
//                   |
//                   |
//                   V



//periodic active-cross model with 1 active particle, a[row][column]
//for 30X30 system we can go upto rho=0.14
//L(50),trelax(4),tmax(2) = 

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          50
# define trelax     10000
# define tmax       120


int main()
{ clock_t tic = clock();

int    a[L][L],a_prime[L][L],position1[L*L]={0},position2[L*L]={0},NN,y,z,i,j,p,q,e;
long int t,tt,Tr,Tm;
double rho,c,d,slot1,slot2,slot3,slot4,block1,block2,rT=1.0,rA,rR,r,dt,pT,pA,pR,pTA,u,v;
char outputfile[200];
init_genrand(87215);
FILE   *fp;


rho=0.05;
rA=15.0;
rR=0.0;



fp = fopen("datafile1.xpm","w");



NN = L*L*rho;
r = 4.0*NN*rT + rA + 2.0*rR;
dt = 1.0/r;
Tr = (int)r*trelax;
Tm = (int)r*tmax;


pT = 4.0*NN*rT/r;
pA = 1.0*rA/r;
pR = 2.0*rR/r;

pTA = pT + pA;

slot1 = 1.0/4.0;
slot2 = 2.0/4.0;
slot3 = 3.0/4.0;
slot4 = 4.0/4.0;




for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0; a_prime[y][z]=0;}
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
        else{a[p][q]=9;  }
              }

      //rotational move of active particle with y+
      else if(a[p][q]==7){
        if(d < 0.5){a[p][q]=5; }
        else{a[p][q]=8; }
              }

     //rotational move of active particle with x-
     else if(a[p][q]==8){
        if(d < 0.5){a[p][q]=7; }
        else{a[p][q]=9; }
              }

     //rotational move of active particle with y-
     else if(a[p][q]==9){
        if(d < 0.5){a[p][q]=5; }
        else{a[p][q]=8; }
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
        else{a[p][q]=9;  }
              }

      //rotational move of active particle with y+
      else if(a[p][q]==7){
        if(d < 0.5){a[p][q]=5; }
        else{a[p][q]=8; }
              }

     //rotational move of active particle with x-
     else if(a[p][q]==8){
        if(d < 0.5){a[p][q]=7; }
        else{a[p][q]=9; }
                   }

     //rotational move of active particle with y-
     else if(a[p][q]==9){
        if(d < 0.5){a[p][q]=5; }
        else{a[p][q]=8; }
                   }

                }//total rotational move ends


for(y=0;y<L;y++)for(z=0;z<L;z++){a_prime[y][z] = 0;}
for(y=0;y<L;y++)for(z=0;z<L;z++){a_prime[y][z] = a[y][z];}

for(y=0;y<L;y++)for(z=0;z<L;z++){
if(a_prime[p][q]==5){a_prime[(p+1)%L][q]=5; a_prime[(p-2+L)%L][q]=5; a_prime[p][(q-1+L)%L]=5; a_prime[p][(q+1)%L]=5;}
else if(a_prime[p][q]==7){a_prime[(p+1)%L][q]=7; a_prime[(p-2+L)%L][q]=7; a_prime[p][(q-1+L)%L]=7; a_prime[p][(q+1)%L]=7;}
else if(a_prime[p][q]==8){a_prime[(p+1)%L][q]=8; a_prime[(p-2+L)%L][q]=8; a_prime[p][(q-1+L)%L]=8; a_prime[p][(q+1)%L]=8;}
else if(a_prime[p][q]==9){a_prime[(p+1)%L][q]=9; a_prime[(p-2+L)%L][q]=9; a_prime[p][(q-1+L)%L]=9; a_prime[p][(q+1)%L]=9;}
                                }





//print------------------------------------------
fprintf(fp,"/* XPM */");            fprintf(fp,"\n");
fprintf(fp,"static char * MyArray[] = {");fprintf(fp,"\n");
fprintf(fp,"\"%d %d 8 1\",",(int)L,(int)L);   fprintf(fp,"\n");      //x,t,type object
fprintf(fp,"\"0      c #ffffff\",");fprintf(fp,"\n");    //for empty
fprintf(fp,"\"2      c #52b7d6\",");fprintf(fp,"\n");    //for passive center
fprintf(fp,"\"3      c #32abd0\",");fprintf(fp,"\n");    //for passive hands
//fprintf(fp,"\"5      c #17202a\",");fprintf(fp,"\n");    //for active center x+ black
fprintf(fp,"\"5      c #fb2232\",");fprintf(fp,"\n");    //for active center x+ red
fprintf(fp,"\"7      c #1bcb0b\",");fprintf(fp,"\n");    //for active center y+ green
fprintf(fp,"\"8      c #f710fa\",");fprintf(fp,"\n");    //for active center x- magenta
fprintf(fp,"\"9      c #f3c72e\",");fprintf(fp,"\n");    //for active center y- yellow
//fprintf(fp,"\"6      c #17202a\",");fprintf(fp,"\n");    //for active hands black
fprintf(fp,"\"6      c #fb2232\",");fprintf(fp,"\n");    //for active hands red

for(y=0;y<L;y++){fprintf(fp,"\""); for(z=0;z<L;z++){fprintf(fp,"%d",a_prime[y][z]); } fprintf(fp,"\","); fprintf(fp,"\n"); } 
//print------------------------------------------


            }//Tm ends


//clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends

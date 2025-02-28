//periodic active-cross model with 1 active particle, a[row][column]
//for 30X30 system we can go upto rho=0.14
//L(100),trelax(10000),tmax(100000),realz(10),(rho=0.13) = 10'

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          50
# define trelax     10000
# define tmax       10000
# define realz      2000


int main()
{ clock_t tic = clock();

int      a[L][L],a_prime[L][L],position1[L*L]={0},position2[L*L]={0},NN,NN2,y,z,i,j,k,count,p,q,e,gap=50000;
long int t,tt,Tr,Tm,x1,x2,rotation=0;
double   rho,rho2,c,d,d1,slot1,slot2,slot3,slot4,block1,block2,rT=1.0,rA,rR,r,dt,pT,pA,pR,pTA,u,v,y1,y2,slope;
double   *s,*xt,*yt,*xtsd,*ytsd,*msd,*msdav;
char     outputfile[200];
init_genrand((unsigned long)time(NULL));
//init_genrand(89751);
FILE     *f1,*fp;



rA=10.0;

//==============
rR=10.0;
rho=0.13;
//==============

rho2=0.16;


//=============================================
sprintf(outputfile, "RRmsd-%d-%.2f-%.2f",(int)rA,rR,rho2);
fp = fopen(outputfile,"w");
//=============================================





//for(rR=1.0; rR<50.0; rR*=1.7){

rotation=0;



NN = L*L*rho;
r = 4.0*NN*rT + rA + 2.0*rR;
dt = 1.0/r;
Tr = (int)r*trelax;
Tm = (int)r*tmax;




s=(double *)calloc(Tm,sizeof(double));
xt=(double *)calloc(Tm,sizeof(double));
yt=(double *)calloc(Tm,sizeof(double));
xtsd=(double *)calloc(Tm,sizeof(double));
ytsd=(double *)calloc(Tm,sizeof(double));
msd=(double *)calloc(Tm,sizeof(double));
msdav=(double *)calloc(Tm,sizeof(double));
for(t=0;t<Tm;t++){s[t]=0; xt[t]=0; yt[t]=0; xtsd[t]=0; ytsd[t]=0; msd[t]=0; msdav[t]=0;}

slope=0;


pT = 4.0*NN*rT/r;
pA = 1.0*rA/r;
pR = 2.0*rR/r;

pTA = pT + pA;

slot1 = 1.0/4.0;
slot2 = 2.0/4.0;
slot3 = 3.0/4.0;
slot4 = 4.0/4.0;



//================================================
//================================================
//================================================
//================================================

for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0;}
//random initialization of passive particles
  i=0; while (i<NN){
  y=genrand_real1()*L; z=genrand_real1()*L;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){
  a[y][z]=2; a[(y+1)%L][z]=3; a[(y-1+L)%L][z]=3; a[y][(z+1)%L]=3; a[y][(z-1+L)%L]=3;
  position1[i] = y; position2[i] = z;
  i++ ; }
                  }
                
NN2 = L*L*rho2;
                  
                  
k=NN; while (k<NN2){

for(y=0;y<L;y++)for(z=0;z<L;z++){
//choosing particle for thermal move
d = genrand_real1();
i = genrand_real1()*k;
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

                          }//total thermal move ends
                    }
   
   
y=genrand_real1()*L; z=genrand_real1()*L;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){
  a[y][z]=2; a[(y+1)%L][z]=3; a[(y-1+L)%L][z]=3; a[y][(z+1)%L]=3; a[y][(z-1+L)%L]=3;
  position1[k] = y; position2[k] = z;
  k++ ; }

     }            

//randomly choose one active among all passive particles
  i=0; while (i<1){
  y = position1[i]; z = position2[i];
  a[y][z]=7; a[(y+1)%L][z]=6; a[(y-1+L)%L][z]=6; a[y][(z+1)%L]=6; a[y][(z-1+L)%L]=6;
  i++ ;
                  }

//================================================
//================================================
//================================================
//================================================



//f1 = fopen("output","w");for(y=0;y<L;y++){for(z=0;z<L;z++){fprintf(f1,"%d",a[y][z]); } fprintf(f1,"\n"); } fprintf(f1,"\n");


//count=0;for(y=0;y<L;y++){for(z=0;z<L;z++){if(a[y][z]==2 || a[y][z]==7)count++; }}printf("%d\n",count);



//realization loop starts
for(e=0;e<realz;e++){

//=================================================================================================
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
//choosing active particle for rotational move
i = 0;
p = position1[i]; q = position2[i];

      //rotational move of active particle with x+
      if(a[p][q]==5){
        if(d < 0.5){a[p][q]=7; }
        else {a[p][q]=9; }
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
//choosing active particle for rotational move
i = 0;
p = position1[i]; q = position2[i];

      //rotational move of active particle with x+
      if(a[p][q]==5){
        if(d < 0.5){a[p][q]=7;  rotation++; }
        else {a[p][q]=9;  rotation++; }
              }

      //rotational move of active particle with y+
      else if(a[p][q]==7){
        if(d < 0.5){a[p][q]=5;  rotation++; }
        else {a[p][q]=8;  rotation++; }
              }

     //rotational move of active particle with x-
     else if(a[p][q]==8){
        if(d < 0.5){a[p][q]=7;  rotation++; }
        else {a[p][q]=9;  rotation++; }
                   }

     //rotational move of active particle with y-
     else if(a[p][q]==9){
        if(d < 0.5){a[p][q]=5;  rotation++; }
        else {a[p][q]=8;  rotation++; }
                   }

                }//total rotational move ends


       xt[t] = u; yt[t] = v;

       xtsd[t] += xt[t]*xt[t];  ytsd[t] += yt[t]*yt[t];


            }//Tm ends


      }//realz ends

  for(t=0;t<Tm;t++){msdav[t] = xtsd[t] + ytsd[t] ;}


//measuring rotational rate in simulation
//fprintf(fp,"%lf  %lf\n",rR,1.0*rotation/realz/tmax/2.0);


//     }//parameter loop ends



//printing MSD for rotating active tracer
for(t=gap;t<Tm;t*=1.15){fprintf(fp,"%lf   %lf\n",1.0*t*dt,msdav[t]/realz*1.0);}


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
free(s); free(xt); free(yt); free(xtsd); free(ytsd); free(msd); free(msdav);
}//main ends

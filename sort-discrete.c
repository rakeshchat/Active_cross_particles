
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L         30
#define tmax      10



int main()
{

int    NN,time,Tm,lg,counter=1,gap,init,i;  
double rho,rT,rA,rR,r;

FILE *fo;

rho=0.06;
rT=1.0;
rA=10.0;
rR=0.0;

NN = rho*L*L;
r = 4.0*NN*rT + rA + 2.0*rR;



//printf("%lf\n",r);

//printf("Total Frames=%d\n",Tm);



fo = fopen("sort-discrete.sh","w");

Tm = 38250;
/*
fprintf(fo,"awk \'");
   for(time=150;time<Tm+1;time+=150){
   fprintf(fo,"NR==%d",time); 
   if(time!=Tm){fprintf(fo," || "); }
}
fprintf(fo,"\' 1no-rot_data.txt > file1");
*/


i=100;
init = 29400;

Tm   = init + 150*15;


fprintf(fo,"awk \'");
   for(time=init;time<=Tm-150;time+=150){
   fprintf(fo,"NR==%d",time); 
   if(time != (init + 150*14) ){fprintf(fo," || "); }
}
//if(init<=5){fprintf(fo,"\' file2 > rAv-r0p0%d",init);}
//else{fprintf(fo,"\' data-act-norot > rAv-r0p0%d",init+1); }

fprintf(fo,"\' 1_data-act-norot > file%d",i);
//fprintf(fo,"\' 2_data-act-rot > rv-rA%d",i);

fclose(fo);
}//main ends

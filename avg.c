//To do averaging with different data files


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M  10000

int main()
{


int    l,k,NN,file = 10;
double *c1,*c2,a,b;
c1=(double *)calloc(M,sizeof(double));
c2=(double *)calloc(M,sizeof(double));
char inputfile1[501],inputfile2[501],outputfile[501];


FILE *fo,*f2,*f3;


fo = fopen("output","w");


f2 = fopen("datalist","r");
k=0; while (fscanf(f2,"%s\n",inputfile2) !=EOF){

          f3 = fopen(inputfile2,"r");
          l=0; while (fscanf(f3,"%lf   %lf\n",&a,&b) !=EOF){ c1[l]=a; c2[l]+=b; l++ ; }
          NN=l;

          fclose(f3);

k++; }//reading inputfile1 ends


for(l=0;l<NN;l++){fprintf(fo,"%lf    %lf\n",c1[l],c2[l]/file*1.0); }



fclose(f2); fclose(fo);
free(c1); free(c2);
}//main ends

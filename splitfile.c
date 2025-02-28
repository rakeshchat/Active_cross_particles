//for spliting files
//this program creates a file with all the splitted file names.


//to generate gif:
//convert -delay 30 -loop 0 *.xpm animated.gif
//convert -delay 100 -loop 1 *.xpm animated.gif
//(20 means 20/100 second delay between each frame, 1 means no loop)
//gifsicle -i input.gif --optimize=3 -o output.gif
//convert input.gif -fuzz 10% -layers Optimize output.gif

//ffmpeg -i input.ogv -vcodec libx264 "output.mp4"



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L         50
#define tmax      120



int main()
{

int    NN,time,Tm,lg,counter=1,gap,step=50;  
double rho,rT=1.0,rA,rR,r;

FILE *fo;

rho=0.05;
rA=15.0;
rR=0.0;


NN = rho*L*L;
r = 4.0*NN*rT + rA + 2.0*rR;

Tm = (int)r*tmax;


gap=11;        //gap = 8 for only non rotation, 11 for rotation


lg = L+gap;


fo = fopen("dosplit.sh","w");
for(time=0;time<Tm;time+=step){
fprintf(fo,"awk \'NR>=%d && NR<=%d\' datafile1.xpm > xfile%04d.xpm\n",lg*time + 1,lg + lg*time,counter); counter++;
                     }


printf("Total Frames= %.1f\n",1.0*Tm/step);


fclose(fo);
}//main ends

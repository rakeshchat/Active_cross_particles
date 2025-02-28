//for spliting files
//this program creates a file with all the splitted file names.


//to generate gif:
//convert -delay 100 -loop 0 *.xpm animated.gif
//convert -delay 20 -loop 1 *.xpm animated.gif
//(20 means 20/100 second delay between each frame, 1 means and no loop)
//gifsicle -i input.gif --optimize=3 -o output.gif
//convert input.gif -fuzz 10% -layers Optimize output.gif

//ffmpeg -i input.ogv -vcodec libx264 "output.mp4"



#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main()
{

int    NN,i;

FILE *fo;


//put this number
NN = 42750;



fo = fopen("sort-continuous.sh","w");


for(i=0;i<15;i++){
fprintf(fo,"awk \'NR>=%d + 1 + 150*%d && NR<=%d + 150 + 150*%d \' 2_data-act-rot > file%d\n",NN,i,NN,i,i);
                 }




fclose(fo);
}//main ends

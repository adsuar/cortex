#include <stdio.h>
#include <string.h>
#include <stdlib.h>

main(int argc,char **argv)
{
   if(argc != 2) exit(-1);
   int max = atoi(argv[1]);

   for(int i=2;i<max;i++)
      if(max%i == 0)
         printf("%d\n",i);
}

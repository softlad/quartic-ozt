#include <stdio.h>
#include <stdlib.h>

/****************************************************************************/
void *MAT_malloc(size_t size)
{
  void *ptr;
if (size!=0)
{
  if ((ptr=malloc(size))==NULL)
  {
	 printf("\nOut Of Memory in MAT_malloc\n");
	// getch();
	 exit(0);
  }
}
else ptr=NULL;
  return(ptr);
}

/****************************************************************************/
void *MAT_realloc(void* ptr,size_t size)
{

  if ((ptr=realloc(ptr,size))==NULL)
  {
	 printf("\nOut Of Memory in MAT_realloc\n");
	 //getch();
	 exit(0);
  }
  return(ptr);
}
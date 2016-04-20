/*Laplace Approximation
 * 
 *
               File : laplace.c
            Authors : Susheel Sagar , Biswajeet Mohanty
 *          Created : 2015-04-10
 *    Last Modified : 2015-05-30
 * Last Modified by : Susheel Sagar
 * 
****************************************************************************/

/* 
 * Compile with:
 * mpicc -o lp laplace.c 
 */
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#define SIZE 1024/* assumption: SIZE a multiple of number of nodes */
                        /* Hint: use small sizes when testing, e.g., SIZE 8 */
#define FROM_MASTER 1   /* setting a message type */
#define FROM_WORKER 2   /* setting a message type */
#define DEBUG 0         /* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
int check=0;
static void print_matrix_A(void)
{
    int i, j;
      printf("\n\n*************printing A************\n\n");
	for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
            printf(" %7.2f", a[i][j]);
        printf("\n\n");
    	}
}

static void print_matrix_B(void)
{
    int i, j;
	 printf("\n\n*************printing B************\n\n");
      for (i = 0; i < SIZE; i++)
    {
        for (j = 0; j < SIZE; j++)
            printf(" %7.2f", b[i][j]);
        printf("\n\n");
    }
}


int validate(double k)
{
	double temp,temp2;
	int i,j;
       	for(i=1;i<SIZE-1;i++)
       	for(j=1;j<SIZE-1;j++)
	{	
//		printf("\ni= %d, j= %d",i,j);
		if(check==0)
		{temp=fabs(a[i][j]-b[i][j]);	}	
//		printf("\n temp when check=0----->%lf",temp);
		else{
		temp = fabs(b[i][j]-a[i][j]);}
//		printf("\n temp when check=1----> %lf",temp);	
	if(temp2<temp)
	temp2=temp;
		
	}
	if(temp2<k){printf("\ntemp2= %lf",temp2);
	return 0;}
	else return 1;
	
}

int main(int argc, char** argv)
{
	int i,j,mtype,dest,src,rows,row_offset,iteration=1,proc_count,core,first=0;
	char filename[40];double accpt_value,temp,start_time,end_time; FILE *fp;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_count);
        MPI_Comm_rank(MPI_COMM_WORLD,&core);
	if(core==0)
        { /*MASTER TASK*/
	if(first==0)
	printf("\n enter the file name to take as input to the matrix"); /* reading the input file path */
        scanf("%s",filename);
        printf("\nenter the acceptance value"); /*reading the acceptance value*/
        scanf("%lf",&accpt_value);
        fp = fopen(filename,"r"); /* reading the file*/
        for(i=0;i<SIZE;i++)      /* loading the data from the file into the matrix a */
        {
                for(j=0;j<SIZE;j++)
                {
                        fscanf(fp,"%lf",&temp);
                        a[i][j] = temp;
                        fscanf(fp,"\t"); /* assuming that each column in a matrix is separated by a tab*/
                }
                fscanf(fp,"\n");  /*assuming that each row in a matrix is separated by a line*/
        } 
        
  //      print_matrix_A();

		printf("\nSIZE = %d, number of nodes = %d\n", SIZE, proc_count);
                start_time = MPI_Wtime();
		
		for(j=0;j<SIZE;j++)
		{ b[0][j] = a[0][j];b[SIZE-1][j]=a[SIZE-1][j];	}
		for(i=0;i<SIZE;i++)
		{b[i][0]=a[i][0];b[i][SIZE-1]=a[i][SIZE-1];	}
         
		if(proc_count==1)
		{
			while(iteration==1)
			{
				if(check==0)
				{ printf("in check = 0");
					for(i=1;i<SIZE-1;i++)
					for(j=1;j<SIZE-1;j++){
					b[i][j]= (a[i-1][j] + a[i+1][j] + a[i][j-1] + a[i][j+1])/4;
					printf("\n");
					}
				}
				else
				{ printf(" in check = 1");
					for(i=1;i<SIZE-1;i++)
                                        for(j=1;j<SIZE-1;j++){
					a[i][j]= (b[i-1][j] + b[i+1][j] + b[i][j-1] + b[i][j+1])/4;		
					printf("\n");
					}
				}
				iteration=validate(accpt_value);
				if(check==0) check=1;
				else check=0;
				 
			}	
		} /*end of code for single node*/
		
		else
		{ 
			while(iteration==1)
			{
			if(check==0)
			{
				printf("\ncheck value is---->%d\n",check);
				rows= SIZE/proc_count;
				mtype = FROM_MASTER;
				row_offset=rows;
				for(dest=1;dest<proc_count;dest++)	
				{	printf("destination == %d",dest);					
					MPI_Send(&row_offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
					MPI_Send(&check, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
    				        MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
					if(dest!=proc_count-1)
					MPI_Send(&a[row_offset-1][0], (rows+2)*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
					else
					MPI_Send(&a[row_offset-1][0],(rows+1)*SIZE, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
					MPI_Send(&b[row_offset][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
					row_offset+=rows;	
				}
				/*Let MASTER do its part of the work*/
				for(i=1;i<rows;i++)
				for(j=1;j<SIZE-1;j++){
				b[i][j]= (a[i-1][j] + a[i+1][j] + a[i][j-1] + a[i][j+1])/4;}
				//printf("\nb[i][j]---->%lf",b[i][j]);
				/*Receiving from the worker*/
				mtype= FROM_WORKER;
				row_offset=rows;
				for(src=1;src<proc_count;src++)
				{
					MPI_Recv(&b[row_offset][0], rows*SIZE, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD,&status);	
					row_offset+=rows;
				}
				
			}
			else if(check==1)
			{
				printf("\ncheck value is -----> %d",check);
				rows= SIZE/proc_count;
				mtype= FROM_MASTER;
				row_offset = rows;
				for(dest=1;dest<proc_count;dest++)
				{	
					printf("\nsending data");
        			        MPI_Send(&row_offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                	                MPI_Send(&check, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				        MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
					if(dest!=proc_count-1)
					MPI_Send(&b[row_offset-1][0],(rows+2)*SIZE,MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
					else
                                        MPI_Send(&b[row_offset-1][0],(rows+1)*SIZE, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
                                        MPI_Send(&a[row_offset][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                                        
                                        row_offset+=rows; 
				}
				/*Let MASTER do its part of the work*/
				for(i=1;i<rows;i++)
                                for(j=1;j<SIZE-1;j++){
                                a[i][j]= (b[i-1][j] + b[i+1][j] + b[i][j-1] + b[i][j+1])/4;}
				//printf("\n\nthe value of a[i][j] is %lf",a[i][j]);	
				 /*Receiving from the worker*/
                                mtype= FROM_WORKER;
                                row_offset=rows;
                                for(src=1;src<proc_count;src++)
                                {
                                        MPI_Recv(&a[row_offset][0], rows*SIZE, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD,&status);
                                        row_offset+=rows;
                                }
			}
			iteration=validate(accpt_value);
			printf("\niteration value== %d",iteration);
			if(check==1) check=0;
			else check=1;
			
		
			}	/*end of while loop*/			
		}	/*end of multiple nodes code*/	
	end_time = MPI_Wtime();
	//if(check==1)
	//print_matrix_A();
//	 print_matrix_B();
	printf("\n\nExecution time on %2d nodes: %f\n", proc_count, end_time-start_time);

	}	/*end of (core==0)  master node code */

	/* entering the worker node code*/
	else
	{
		while(iteration==1)
		{
//		printf("core no is printing ---->>>>>>>>>>>> %d",core);
		mtype = FROM_MASTER;
                MPI_Recv(&row_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD,&status);
                MPI_Recv(&check, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD,&status);
		MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD,&status);
//		printf("check value received is %d",check);
		if(check==0)
		{
			if(core!=proc_count-1)
				MPI_Recv(&a[row_offset-1][0], (rows+2)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
			else
				MPI_Recv(&a[row_offset-1][0], (rows+1)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
			MPI_Recv(&b[row_offset][0], (rows)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
			for(i=row_offset;i<rows+row_offset;i++)
			for(j=1;j<SIZE-1;j++)
			{
				b[i][j] = (a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1])/4;
			//	printf("");
			}
			mtype=FROM_WORKER;
			MPI_Send(&b[row_offset][0],rows*SIZE,MPI_DOUBLE,0,mtype,MPI_COMM_WORLD);
		}

		else if (check==1)
                {
		//	printf("core no is ---> %d",core);
                	if(core!=proc_count-1)
		           { printf("");
			    MPI_Recv(&b[row_offset-1][0], (rows+2)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
				}
                	else
		                MPI_Recv(&b[row_offset-1][0], (rows+1)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
                	MPI_Recv(&a[row_offset][0], (rows)*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&status);
	                for(i=row_offset;i<rows+row_offset;i++)
        	        for(j=1;j<SIZE-1;j++)
                	{
                        	a[i][j] = (b[i-1][j] + b[i+1][j] + b[i][j+1] + b[i][j-1])/4;
			//	printf("xyzzzzzz");
                	}
               		 mtype=FROM_WORKER;
                	MPI_Send(&a[row_offset][0],rows*SIZE,MPI_DOUBLE,0,mtype,MPI_COMM_WORLD);
                }
		}
	}	/* end of worker tasks*/
	
	MPI_Finalize();
	return 0;
}	/* END OF MAIN*/




/*checkerboard version of  Matrix-Matrix multiplication
 * 
 *             File : checkerboard.c
 *        Author(s) : Susheel Sagar, Biswajeet Mohanty
 *          Created : 2015-04-26
 *    Last Modified : 2015-05-22
 * Last Modified by : Susheel Sagar
 * 
 ***************************************************************************/

/* 
 * Compile with:
 * mpicc -o cb checkerboard.c 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 256  /* assumption: SIZE a multiple of number of nodes */
#define FROM_MASTER 1   /* setting a message type */
#define FROM_WORKER 2   /* setting a message type */
#define DEBUG 0         /* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

static void init_matrix(void)
{
    int i,j;
     for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
            /* Simple initialization, which enables us to easily check
             * the correct answer. Each element in c will have the same 
             * value as SIZE after the matmul operation.
             */
            a[i][j] = 2.0;
            b[i][j] = 1.0;
	    c[i][j] = 0.0;
        }
}
static void print_matrix(void)
{
    int i, j;
      for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
            printf(" %7.2f", c[i][j]);
        printf("\n");
    }
}

int main(int argc, char **argv)
{
 int core,proc_count,i,j,k,mtype,row_offset,column_offset,rows,columns,src,dest;
 double start_time, end_time;
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &core);
	if(core==0)
	{
		printf("SIZE = %d, number of nodes = %d\n", SIZE, proc_count); //printing the number of nodes
	        init_matrix();  // initialization of the matrix
        	start_time = MPI_Wtime(); //starting the clock
		
		if(proc_count==1) /* code to run on 1 node*/
		{
			for (i = 0; i < SIZE; i++) {
        		for (j = 0; j < SIZE; j++) {
           		for (k = 0; k < SIZE; k++)
               		 c[i][j] = c[i][j] + a[i][k] * b[k][j];
			}

			}
		}
		else if(proc_count==2) /*code to run on two nodes*/
		{
			mtype = FROM_MASTER;
			rows=SIZE/2;
			row_offset = rows;
			MPI_Send(&rows, 1, MPI_INT, 1 , mtype, MPI_COMM_WORLD);
			MPI_Send(&row_offset, 1, MPI_INT, 1 , mtype, MPI_COMM_WORLD);
			MPI_Send(&a[row_offset][0], (SIZE*SIZE)/2, MPI_DOUBLE,1,mtype,MPI_COMM_WORLD);
			MPI_Send(&b, SIZE*SIZE, MPI_DOUBLE, 1, mtype, MPI_COMM_WORLD);
			
			for (i = 0; i < row_offset; i++) 
			{
		            for (j = 0; j < SIZE; j++) 
				{
                                for (k = 0; k < SIZE; k++)
                		    c[i][j] = c[i][j] + a[i][k] * b[k][j];
            			}
        		}
			 mtype = FROM_WORKER;

			 MPI_Recv(&row_offset, 1, MPI_INT, 1, mtype, MPI_COMM_WORLD, &status);
			  MPI_Recv(&rows, 1, MPI_INT, 1, mtype, MPI_COMM_WORLD, &status);
			 MPI_Recv(&c[row_offset][0], SIZE*SIZE/2, MPI_DOUBLE, 1, mtype, MPI_COMM_WORLD, &status);
							

		}
		else if(proc_count==4)  /*code to run on 4 nodes*/
		{
			mtype= FROM_MASTER;
			rows=SIZE/2; columns=SIZE/2;
			row_offset = 0; column_offset = columns;
			for(dest=1;dest<proc_count;dest++)
			{
				MPI_Send(&rows, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
				MPI_Send(&columns, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
				MPI_Send(&row_offset, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
				MPI_Send(&column_offset, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);	
				MPI_Send(&a[row_offset][0], (SIZE*SIZE)/2, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
				for(i=0;i<SIZE;i++)
				MPI_Send(&b[i][column_offset], SIZE/2, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
				if(dest==1){row_offset+=rows;column_offset=0;}
				if(dest==2){column_offset=columns;}
			}
			for(i=0;i<rows;i++)
			for(j=0;j<columns;j++){
			for(k=0;k<SIZE;k++)
			 c[i][j] = c[i][j] + a[i][k] * b[k][j];
			
			}
			mtype = FROM_WORKER;
			for(src=1;src<4;src++)
			{
				if(src==1){row_offset=0;column_offset=columns;}			
			        if(src==2){row_offset=rows;column_offset=0;}
				if(src==3){row_offset=rows;column_offset=columns;}
			for(i=row_offset;i<row_offset+rows;i++)
			for(j=column_offset;j<column_offset+columns;j++)
			MPI_Recv(&c[i][j], 1, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
			}
		}
		else if(proc_count==8) /*code to run on 8 nodes*/
		{
			mtype = FROM_MASTER;
			rows=SIZE/4; columns=SIZE/2;		
			row_offset=0; column_offset=columns; 
			for(dest=1;dest<proc_count;dest++)
			{
				MPI_Send(&rows, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
                                MPI_Send(&columns, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
                                MPI_Send(&row_offset, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
                                MPI_Send(&column_offset, 1, MPI_INT, dest , mtype, MPI_COMM_WORLD);
                                MPI_Send(&a[row_offset][0], (SIZE*SIZE)/4, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
                                for(i=0;i<SIZE;i++)
                                MPI_Send(&b[i][column_offset], SIZE/2, MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
				if(dest==1){row_offset+=rows;column_offset=0;}
				if(dest==2){column_offset+=columns;}
				if(dest==3){row_offset+=rows;column_offset=0;}
				if(dest==4){column_offset+=columns;}
				if(dest==5){row_offset+=rows;column_offset=0;}
				if(dest==6){column_offset=columns;}
			}
				for(i=0;i<rows;i++)
	                        for(j=0;j<columns;j++){
        	                for(k=0;k<SIZE;k++)
                	        c[i][j] = c[i][j] + a[i][k] * b[k][j];
				}
			mtype = FROM_WORKER;
                        for(src=1;src<8;src++)
                        {
				 if(src==1){row_offset=0;column_offset=columns;}	
				 if(src==2){row_offset+=rows;column_offset=0;}
				 if(src==3){column_offset=columns;}	
				if(src==4){row_offset+=rows;column_offset=0;}
				if(src==5){column_offset=columns;}
				if(src==6){row_offset+=rows;column_offset=0;}
				if(src==7){column_offset=columns;}
				for(i=row_offset;i<row_offset+rows;i++)
                        for(j=column_offset;j<column_offset+columns;j++)
                        MPI_Recv(&c[i][j], 1, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);

			}	
		}
	 end_time = MPI_Wtime();  /* ending the time*/
	 print_matrix();  /* printing the resultant  matrix*/
         printf("\n\nExecution time on %2d nodes: %f\n", proc_count, end_time-start_time); /* printing the execution time */

	}
	else
	{/* worker node code*/
		if(proc_count==2) /* code for 2 nodes*/
		{

			mtype = FROM_MASTER;
			MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&row_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&a[row_offset][0],(SIZE*SIZE)/2 , MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
		        MPI_Recv(&b, SIZE*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
			for (i=row_offset; i<row_offset+rows; i++)
		            for (j=0; j<SIZE; j++) {
			    	for (k=0; k<SIZE; k++)
                    		c[i][j] = c[i][j] + a[i][k] * b[k][j];							
				}
			
			mtype = FROM_WORKER;			
	        	MPI_Send(&row_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
		        MPI_Send(&c[row_offset][0], (SIZE*SIZE)/2, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
			
		}
		if(proc_count==4) /* code for 4 nodes*/
		{
			mtype= FROM_MASTER;
			MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&row_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&column_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&a[row_offset][0],(SIZE*SIZE)/2 , MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
			for(i=0;i<SIZE;i++)
			MPI_Recv(&b[i][column_offset], SIZE/2, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
			mtype = FROM_WORKER;
			for (i=row_offset; i<row_offset+rows; i++)
                            for (j=column_offset; j<column_offset+columns; j++) {
                                for (k=0; k<SIZE; k++)
                                c[i][j] = c[i][j] + a[i][k] * b[k][j];
				MPI_Send(&c[i][j],1,MPI_DOUBLE,0,mtype,MPI_COMM_WORLD);
                                }
			
			
			
			
		}
		if(proc_count==8)  /*code for 8 nodes*/
		{
			 mtype= FROM_MASTER;
                        MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&row_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&column_offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                        MPI_Recv(&a[row_offset][0],(SIZE*SIZE)/4 , MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
                        for(i=0;i<SIZE;i++)
                        MPI_Recv(&b[i][column_offset], SIZE/2, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
                        mtype = FROM_WORKER;
                        for (i=row_offset; i<row_offset+rows; i++)
                            for (j=column_offset; j<column_offset+columns; j++) {
                                for (k=0; k<SIZE; k++)
                                c[i][j] = c[i][j] + a[i][k] * b[k][j];
                                MPI_Send(&c[i][j],1,MPI_DOUBLE,0,mtype,MPI_COMM_WORLD);
                                }

		}
	}
	 MPI_Finalize();
         return 0;

}




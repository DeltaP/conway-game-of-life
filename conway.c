/*
 * Gregory Petropoulos
 *
 * This is my ping pong program for assignment 4
 *
 * To compile:  mpicc -g -Wall -std=c99 -o CGL conway.c
 * To run:  mpiexec -n 1 ./CGL <input file name> <partition> <iteration> <interval>
 *          <> -> mandatory
 *          [] -> optional
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>

// -----------------------------------------------------------------
// global variables
int nprocs;
int my_rank;
// Local physical field size
int field_width;        // Width and height of field on this processor
int field_height;       // (should be local_width+2, local_height+2)
int local_width;
int local_height;
int width;
int height;
int ncols;
int nrows;
int *field_a;      // The local data fields
int *field_b;
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Ends the program on an error and prints message to node 0
void cleanup (int my_rank, const char *message) {
  if (my_rank==0) {printf("%s\n",message);}
  MPI_Finalize();                                       /* kills mpi                            */
  exit(0);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// reads in file
bool fileread (char *filename, char *partition) {

  // Open the file
  if( my_rank==0 )
    printf( "Opening file %s\n", filename );
  FILE *fp = fopen( filename, "r" );
  if( !fp ) {
    printf( "Error: The file '%s' could not be opened.\n", filename );
    return false;
  }

  // Read the PGM header, which looks like this:
  //  |P5                magic version number
  //  |900 900           width height
  //  |255               depth
  char header[10];
  int  depth;
  int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
  if( rv != 4 ) {
    if(my_rank==0) 
      printf( "Error: The file '%s' did not have a valid PGM header\n", filename );
      return false;
  }
  if( my_rank==0 ){
    printf( "%s: %s %i %i %i\n", filename, header, width, height, depth );
  }
  // Make sure the header is valid
  if( strcmp( header, "P5") ) {
    if(my_rank==0) 
      printf( "Error: PGM file is not a valid P5 pixmap.\n" );
    return false;
  }
  if( depth != 255 )
  {
    if(my_rank==0) 
      printf( "Error: PGM file has depth=%i, require depth=255 \n", depth );
    return false;
  }

  //
  // Make sure that the width and height are divisible by the number of
  // processors in x and y directions
  //
  if (strcmp(partition, "slice") == 0) {
    ncols = 1;
    nrows = nprocs;
  }
  else if (strcmp(partition, "checkerboard")==0) {
    ncols = sqrt(nprocs);
    nrows = sqrt(nprocs);
  }
  else {
    // Create the array!
    local_width  = width;
    local_height = height;
    field_width  = width  + 2;
    field_height = height + 2;
    field_a = (int *)malloc( field_width * field_height * sizeof(int));
    field_b = (int *)malloc( field_width * field_height * sizeof(int));
    int b; 
    int ll;
    for (int y = 0; y < field_height; y++ ) {
      for(int x = 0; x < field_width; x++ ) {
        if ((x == 0) || (y == 0) || (x == field_width - 1) || (y == field_height - 1)) {
          ll = (y * field_width + x);
          field_a[ll] = 0;
          field_b[ll] = 0;
        }
        else {
          // Read the next character
          b = (int)fgetc( fp );
          if( b == EOF ) {
            printf( "Error: Encountered EOF at [%i,%i]\n", y,x );
            return false;
          }
           
          // From the PGM, black cells (b=0) are bugs, all other 
          // cells are background 
          b = (b==0)?1:0;

          ll = (y * field_width + x);
          field_a[ ll ] = b;
          field_b[ ll ] = b;
        }
      }
    }

    fclose(fp);
    return true;
  }

  fclose(fp);

  if (width % ncols) {
    if (my_rank==0)
      printf( "Error: %i pixel width cannot be divided into %i cols\n", width, ncols );
    return false;
  }
  if (height % nrows) {
    if (my_rank==0)
      printf( "Error: %i pixel height cannot be divided into %i rows\n", height, nrows );
    return false;
  }

  // Divide the total image among the local processors
  local_width = width / ncols;
  local_height = height / nrows;

  // Create the array!
  field_width = local_width + 2;
  field_height = local_height + 2;
  field_a = (int *)malloc( field_width * field_height * sizeof(int));
  field_b = (int *)malloc( field_width * field_height * sizeof(int));
  char *temp=(char *)malloc( local_width * local_height * sizeof(char));


	int i, j;
	MPI_File	mpi_infile;
  MPI_Status status;
	int rc;
  //int temp_size=local_row*local_col;
  //char *temp=(char*)malloc(temp_size*sizeof(char));

  MPI_Aint extent;
  MPI_Datatype etype, filetype, contig;
  
  etype = MPI_CHAR;
  MPI_Type_contiguous(local_width, etype, &contig);

  MPI_Type_contiguous(2, etype, &contig); 
  extent = 6 * sizeof(int); 
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype); 

  return true;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// counts the number of bugs
void measure (int iteration) {
  int local_sum = 0;
  int global_sum = 0;

  for (int j = 0; j < local_height; j++) {
    for (int i = 0; i < local_width; i++) {
      local_sum += (iteration%2==0) ? (field_a[(j+1)*(field_width)+(i+1)]) : (field_b[(j+1)*(field_width)+(i+1)]);
    }
  }
  
  if (nprocs > 1) {MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);}
  else {global_sum = local_sum;}

  if (my_rank == 0) {printf("BUGCOUNT for iteration %i = %i \n", iteration, global_sum);}
}
// -----------------------------------------------------------------


  
// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  char in_file[100];
  char partition[100];
  int iterations, interval, neighbor;

  MPI_Init(&argc, &argv);                               /* start up MPI                         */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);               /* get the number of processes          */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes  */

  if (argc < 5) {
    cleanup(my_rank, "Error:  Too few arguments");
  }
  else if (argc == 5) {
    strcpy(in_file,   argv[1]);
    strcpy(partition, argv[2]);
    if (!((strcmp(partition, "slice") == 0) || (strcmp(partition, "checkerboard") == 0) || (strcmp(partition, "none") == 0))) {
      cleanup(my_rank, "Error:  Incorrect partition option.  Enter 'slice', 'checkerboard', or 'none'");
    }
    iterations = atoi(argv[3]);
    interval   = atoi(argv[4]);
  }
  else if (argc > 5) {
    cleanup(my_rank, "Error:  Too many arguments");
  }

  bool go_on = fileread(in_file, partition);
  if (go_on == false) {cleanup(my_rank, "Error:  fileread returned false, quitting program");}

  for (int i = 0; i < iterations; i++) {
    if (i%interval == 0) {measure(i);}
    for (int y = 0; y < local_height; y++) {
      for (int x = 0; x < local_width; x++) {
        int yb = y+1;
        int xb = x+1; 
        neighbor = 0;
        if (i%2 == 0) {
          neighbor += field_a[(yb-1)*field_width+xb+1];
          neighbor += field_a[(yb-1)*field_width+xb];
          neighbor += field_a[(yb-1)*field_width+xb-1];
          neighbor += field_a[(yb)*field_width+xb+1];
          neighbor += field_a[(yb)*field_width+xb-1];
          neighbor += field_a[(yb+1)*field_width+xb+1];
          neighbor += field_a[(yb+1)*field_width+xb];
          neighbor += field_a[(yb+1)*field_width+xb-1];

          field_b[yb*field_width+xb]=field_a[yb*field_width+xb];

          if ((((neighbor < 2) || (neighbor > 3)) && (field_a[yb*field_width+xb])) || ((neighbor == 3) && (!field_a[yb*field_width+xb]))) {
            field_b[yb*field_width+xb] ^= 1;
          }
        }
        else {
          neighbor += field_b[(yb-1)*field_width+xb+1];
          neighbor += field_b[(yb-1)*field_width+xb];
          neighbor += field_b[(yb-1)*field_width+xb-1];
          neighbor += field_b[(yb)*field_width+xb+1];
          neighbor += field_b[(yb)*field_width+xb-1];
          neighbor += field_b[(yb+1)*field_width+xb+1];
          neighbor += field_b[(yb+1)*field_width+xb];
          neighbor += field_b[(yb+1)*field_width+xb-1];

          field_a[yb*field_width+xb]=field_b[yb*field_width+xb];

				  if ((((neighbor < 2) || (neighbor > 3)) && (field_b[yb*field_width+xb])) || ((neighbor == 3) && (!field_b[yb*field_width+xb]))) {
            field_a[yb*field_width+xb] ^= 1;
          }
        }
      }     // x
    }       // y
  }         // iteration
  cleanup(my_rank, "Thank you for playing!");
  return 0;
}

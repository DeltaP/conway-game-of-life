/*
 * Gregory Petropoulos
 *
 * This is my Conway's Game of Life program for the midterm
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
// global variables --> perhaps make some of these local?
int nprocs;
int my_rank;
int field_width;                                        /* local board size                     */
int field_height; 
int local_width;                                        /* local size of data                   */
int local_height;
int width;                                              /* The total dimension of the field     */
int height;
int ncols;
int nrows;
int *field_a;                                           /* The local data fields                */
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
void fileread (char *filename, char *partition) {
  if (strcmp(partition, "slice") == 0) {                /* determines the data decomposition    */
    ncols = 1;
    nrows = nprocs;
  }
  else if (strcmp(partition, "checkerboard")==0) {
    ncols = sqrt(nprocs);
    nrows = sqrt(nprocs);
  }
  else {
    ncols = 1;
    nrows = 1;
  }

  MPI_File fh;                                          /* sets up MPI input                  */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  int offset =0;
  int newline_count = 0;
  char header[100];
  char read;

  do {                                                  /* reads in the header                */
    MPI_File_read_at_all(fh, offset, &read, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    header[offset]=read;
    if (read == '\n') newline_count++;
    offset++;
    if (offset == 100) cleanup(my_rank, "Error:  Header exceeds 100 characters, check file or recompile code");
  } while (newline_count < 3);
    
  char head[10];                                        /* parses the header and error      */
  int depth;                                            /*   checks                         */
  int rv =                sscanf(header, "%6s\n%i %i\n%i\n", head, &width, &height, &depth);
  if (rv != 4)            cleanup(my_rank,"Error: The file did not have a valid PGM header");
  if (my_rank==0)         printf( "%s: %s %i %i %i\n", filename, head, width, height, depth );
  if (strcmp(head, "P5")) cleanup(my_rank,"Error: PGM file is not a valid P5 pixmap");
  if( depth != 255 )      cleanup(my_rank,"Error: PGM file requires depth=255");
  if (width % ncols)      cleanup(my_rank,"Error: pixel width cannot be divided evenly into cols");
  if (height % nrows)     cleanup(my_rank,"Error: %i pixel height cannot be divided into %i rows\n");

  local_width = width / ncols;                          /* determines the size of           */
  local_height = height / nrows;                        /* the local data                   */

  field_width = local_width + 2;                        /* creates an array with room for   */
  field_height = local_height + 2;                      /* ghosts and boarders              */
  // playing fields
  field_a = (int *)malloc(field_width * field_height * sizeof(int));
  field_b = (int *)malloc(field_width * field_height * sizeof(int));
  // array for the file read
  char *temp=(char *)malloc( local_width * local_height * sizeof(char));

  MPI_Aint extent;                                      /* declares the extent              */
  MPI_Datatype etype, filetype, contig;                 /* derrived data types for IO       */
  MPI_Offset disp = offset;                             /* the initial displacement of      */
                                                        /*   the header                     */
  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;
  
  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR       */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                  */
  extent = width * sizeof(char);                        /* total size of repeatable unit    */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data */
  // reads in the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);

  int b, ll; 
  for (int y = 0; y < field_height; y++ ) {             /* loops through the field          */
    for(int x = 0; x < field_width; x++ ) {
      ll = (y * field_width + x);                       /* puts zeros at the boarders       */
      if ((x == 0) || (y == 0) || (x == field_width - 1) || (y == field_height - 1)) {
        field_a[ll] = 0;
        field_b[ll] = 0;
      }
      else {
        b = (int)temp[x-1+(y-1)*local_width];           /* finds data from read             */
        b = (b==0)?1:0;                                 /* black = bugs; other = no bug     */
        field_a[ll] = b;
        field_b[ll] = b;
      }
    }
  }
  MPI_File_close(&fh);
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
// swaps ghost rows
void summonspectre(int iteration) {
  MPI_Datatype col;                                     /* makes col data type              */
  MPI_Type_vector(field_height, 1, field_width, MPI_INT, &col);
  MPI_Type_commit(&col);

  MPI_Datatype row;                                     /* makes row data type              */
  MPI_Type_vector(field_width, 1, 1, MPI_INT, &row);
  MPI_Type_commit(&row);

  int *pointer = (iteration%2==0)?field_a:field_b;      /* switches between field_a and     */
                                                        /*   field_b                        */
  if ((my_rank%ncols)%2 == 1) {                         /* send right odd                   */
    MPI_Recv(pointer, 1, col, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank%ncols != ncols -1) {             
    MPI_Send((local_width+pointer), 1, col, my_rank + 1, 0, MPI_COMM_WORLD);
  }

  if ((my_rank%ncols)%2 == 1) {                         /* send left odd                    */
    MPI_Send((1+pointer), 1, col, my_rank-1, 0, MPI_COMM_WORLD);
  }
  else if (my_rank%ncols != ncols -1) {
    MPI_Recv((local_width+1+pointer), 1, col, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (my_rank%ncols != 0 && (my_rank%ncols)%2 == 0) {   /* send right even                  */
    MPI_Recv(pointer, 1, col, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank%ncols != ncols -1 && (my_rank%ncols)%2==1) {
    MPI_Send((local_width+pointer), 1, col, my_rank + 1, 0, MPI_COMM_WORLD);
  }

  if (my_rank%ncols != 0 && (my_rank%ncols)%2 == 0) {   /* send left even                   */
    MPI_Send((1+pointer), 1, col, my_rank-1, 0, MPI_COMM_WORLD);
  }
  else if (my_rank%ncols != ncols -1 && (my_rank%ncols)%2==1) {
    MPI_Recv((local_width+1+pointer), 1, col, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if ((my_rank/ncols)%2 == 1) {                         /* send down odd                    */
    MPI_Recv(pointer, 1, row, my_rank - ncols, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank/ncols != nrows -1) {
    MPI_Send((field_width*local_height+pointer), 1, row, my_rank + ncols, 0, MPI_COMM_WORLD);
  }

  if ((my_rank/ncols)%2 == 1) {                         /* send up odd                      */
    MPI_Send((field_width+pointer), 1, row, my_rank-ncols, 0, MPI_COMM_WORLD);
  }
  else if (my_rank/ncols != nrows -1) {
    MPI_Recv((field_width*(local_height+1)+pointer), 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (my_rank/ncols != 0 && (my_rank/ncols)%2 == 0) {   /* send down even                   */
    MPI_Recv(pointer, 1, row, my_rank - ncols, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank/ncols != nrows -1 && (my_rank/ncols)%2==1) {
    MPI_Send((field_width*local_height+pointer), 1, row, my_rank + ncols, 0, MPI_COMM_WORLD);
  }

  if (my_rank/ncols != 0 && (my_rank/ncols)%2 == 0) {   /* send up even                     */
    MPI_Send((field_width+pointer), 1, row, my_rank-ncols, 0, MPI_COMM_WORLD);
  }
  else if (my_rank/ncols != nrows -1 && (my_rank/ncols)%2==1) {
    MPI_Recv((field_width*(local_height+1)+pointer), 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}
// -----------------------------------------------------------------


  
// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  char in_file[1000];
  char partition[100];
  int iterations, interval, neighbor;

  MPI_Init(&argc, &argv);                               /* start up MPI                         */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);               /* get the number of processes          */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes  */

  if (argc < 5) {                                       /* parses command line arguments        */
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

  fileread(in_file, partition);                         /* function call to read the file in    */

  for (int i = 0; i < iterations; i++) {
    summonspectre(i);                                   /* exchanges ghost fields               */
    if (i%interval == 0) {measure(i);}
    for (int y = 0; y < local_height; y++) {            /* loops through the local data         */
      for (int x = 0; x < local_width; x++) {
        int yb = y+1;                                   /* shifts needed becuase the            */
        int xb = x+1;                                   /*   is bigger due to ghost rows        */
        neighbor = 0;                                   /* initialize the neighbor count        */
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

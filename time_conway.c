/*
 * Gregory Petropoulos
 *
 * This is my Conway's Game of Life program for the midterm
 * This is the copy for timing tests
 *
 * To compile:  mpicc -g -Wall -std=c99 -o CGL time_conway.c
 * To run:  mpiexec -n 1 ./CGL <input file name> <partition> <iteration> <m_interval> <w_interval> <t_interval>
 *          <> -> mandatory
 *          [] -> optional
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "mpi.h"

// -----------------------------------------------------------------
// global variables --> perhaps make some of these local?
int nprocs;
int my_rank;
int field_width;                                        /* local board size                       */
int field_height; 
int local_width;                                        /* local size of data                     */
int local_height;
int width;                                              /* The total dimension of the field       */
int height;
int ncols;
int nrows;
int *field_a;                                           /* The local data fields                  */
int *field_b;
char header[100];
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Ends the program on an error and prints message to node 0
void cleanup (int my_rank, const char *message) {
  if (my_rank == 0) printf("%s\n",message);
  MPI_Finalize();                                       /* kills mpi                              */
  exit(0);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// reads in file
void fileread (char *filename, char *partition, int *offset) {
  if (strcmp(partition, "slice") == 0) {                /* determines the data decomposition      */
    ncols = 1;
    nrows = nprocs;
  }
  else if (strcmp(partition, "checkerboard") == 0) {
    ncols = sqrt(nprocs);
    nrows = sqrt(nprocs);
  }
  else {
    ncols = 1;
    nrows = 1;
  }

  MPI_File fh;                                          /* sets up MPI input                      */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  int newline_count = 0;
  char read;
  *offset = 0;

  do {                                                  /* reads in the header                    */
    MPI_File_read_at_all(fh, *offset, &read, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    header[*offset]=read;
    if (read == '\n') newline_count++;
    (*offset)++;
    if (*offset == 100) cleanup(my_rank, "Error:  Header exceeds 100 characters, check file or recompile code");
  } while (newline_count < 3);
    
  char head[10];                                        /* parses the header and error            */
  int depth;                                            /*   checks                               */
  int rv               =  sscanf(header, "%6s\n%i %i\n%i\n", head, &width, &height, &depth);
  if (rv != 4)            cleanup(my_rank,"Error: The file did not have a valid PGM header");
  if (my_rank==0)         printf( "%s: %s %i %i %i\n", filename, head, width, height, depth );
  if (strcmp(head, "P5")) cleanup(my_rank,"Error: PGM file is not a valid P5 pixmap");
  if( depth != 255 )      cleanup(my_rank,"Error: PGM file requires depth=255");
  if (width % ncols)      cleanup(my_rank,"Error: pixel width cannot be divided evenly into cols");
  if (height % nrows)     cleanup(my_rank,"Error: %i pixel height cannot be divided into %i rows\n");

  local_width = width / ncols;                          /* determines the size of                 */
  local_height = height / nrows;                        /* the local data                         */
  field_width = local_width + 2;                        /* creates an array with room for         */
  field_height = local_height + 2;                      /* ghosts and boarders                    */
  // playing fields
  field_a = (int *)malloc(field_width * field_height * sizeof(int));
  field_b = (int *)malloc(field_width * field_height * sizeof(int));
  // array for the file read
  char *temp=(char *)malloc( local_width * local_height * sizeof(char));
  MPI_Aint extent;                                      /* declares the extent                    */
  MPI_Datatype etype, filetype, contig;                 /* derrived data types for IO             */
  MPI_Offset disp = *offset;                            /* the initial displacement of            */
                                                        /*   the header                           */
  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;
  
  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR             */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                        */
  extent = width * sizeof(char);                        /* total size of repeatable unit          */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data       */
  // reads in the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);

  int x, y, b, ll; 
  for (y = 0; y < field_height; y++ ) {                 /* loops through the field                */
    for(x = 0; x < field_width; x++ ) {
      ll = (y * field_width + x);                       /* puts zeros at the boarders             */
      if ((x == 0) || (y == 0) || (x == field_width - 1) || (y == field_height - 1)) {
        field_a[ll] = 0;
        field_b[ll] = 0;
      }
      else {
        b = (int)temp[x-1+(y-1)*local_width];           /* finds data from read                   */
        b = (b==0)?1:0;                                 /* black = bugs; other = no bug           */
        field_a[ll] = b;
        field_b[ll] = b;
      }
    }
  }
  MPI_File_close(&fh);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// writes the game board out to file
void filewrite (char *in_file, int iteration, int offset) {
  char filename[1000];
  sprintf(filename, "%d_", iteration); 
  strcat(filename, in_file);

  MPI_File fh;                                          /* sets up MPI input                      */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  if (my_rank == 0) {                                   /* writes the header                      */
    MPI_File_write(fh,header,offset,MPI_CHAR,MPI_STATUS_IGNORE); 
  }

  char *temp=(char *)malloc( local_width * local_height * sizeof(char));
  int *field_pointer = (iteration%2==0) ? (field_a) : (field_b);

  int x, y, b, ll;
  for (y = 0; y < field_height; y++ ) {                 /* loops through the field                */
    for(x = 0; x < field_width; x++ ) {
      if ((x != 0) && (y != 0) && (x != field_width - 1) && (y != field_height - 1)) {
        ll = (y * field_width + x);                     /* puts zeros at the boarders             */
        b = field_pointer[ll];
        b = (b==0)?0xFF:0;                              /* black = bugs; other = no bug           */
        temp[x-1+(y-1)*local_width] = (char)b;
      }
    }
  }

  MPI_Aint extent;                                      /* declares the extent                    */
  MPI_Datatype etype, filetype, contig;                 /* derrived data types for IO             */
  MPI_Offset disp = offset;                             /* the initial displacement of            */

  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;

  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR             */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                        */
  extent = width * sizeof(char);                        /* total size of repeatable unit          */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data       */
  // writes the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// counts the number of bugs
void measure (int iteration) {
  int local_sum = 0;
  int global_sum = 0;
  int i, j;

  int *pointer = (iteration%2 == 0) ? field_a : field_b;

  for (j = 0; j < local_height; j++) {
    for (i = 0; i < local_width; i++) {
      local_sum += pointer[(j+1)*(field_width)+(i+1)];
    }
  }
  
  if (nprocs > 1) MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  else global_sum = local_sum;

  if (my_rank == 0) printf("BUGCOUNT for iteration %i = %i \n", iteration, global_sum);
}
// -----------------------------------------------------------------
// jijioasdlfk   couch!  for doing that.f

// -----------------------------------------------------------------
// swaps ghost rows
void summonspectre(int iteration, MPI_Request* send, MPI_Request* recv) {
  MPI_Datatype col;                                     /* makes col data type                    */
  MPI_Type_vector(field_height, 1, field_width, MPI_INT, &col);
  MPI_Type_commit(&col);

  MPI_Datatype row;                                     /* makes row data type                    */
  MPI_Type_vector(field_width, 1, 1, MPI_INT, &row);
  MPI_Type_commit(&row);

  /* POINTER LOCATIONS *
   *                   *
   *  |a|b| ... |c|d|  *
   *  |e|              *
   *  |.|              *
   *  |.|              *
   *  |.|              *
   *  |f|              *
   *  |g|              *
   *                   *
   *********************/
  int *pointer = (iteration%2==0)?field_a:field_b;      /* switches between field_a and           */
                                                        /*   field_b                              */
  int *p_a = pointer;                                   /* pointers for communication             */
  int *p_b = pointer + 1;
  int *p_c = pointer + local_width;
  int *p_d = pointer + local_width + 1;
  int *p_e = pointer + field_width;
  int *p_f = pointer + field_width * local_height;
  int *p_g = pointer + field_width * (local_height + 1);



  if ((my_rank%ncols)%2 == 1) {                         /* send left-right                        */
    MPI_Irecv(p_a, 1, col, my_rank-1, 0, MPI_COMM_WORLD, &recv[0]);
    MPI_Isend(p_b, 1, col, my_rank-1, 0, MPI_COMM_WORLD, &send[0]);
    if (my_rank%ncols != ncols -1) {
      MPI_Isend(p_c, 1, col, my_rank+1, 0, MPI_COMM_WORLD, &send[1]);
      MPI_Irecv(p_d, 1, col, my_rank+1, 0, MPI_COMM_WORLD, &recv[1]);
    }
  }
  else if (my_rank%ncols != ncols -1) {             
    MPI_Isend(p_c, 1, col, my_rank+1, 0, MPI_COMM_WORLD, &send[2]);
    MPI_Irecv(p_d, 1, col, my_rank+1, 0, MPI_COMM_WORLD, &recv[2]);
    if (my_rank%ncols != 0) {
      MPI_Irecv(p_a, 1, col, my_rank-1, 0, MPI_COMM_WORLD, &recv[3]);
      MPI_Isend(p_b, 1, col, my_rank-1, 0, MPI_COMM_WORLD, &send[3]);
    }
  }
  MPI_Waitall(8, recv, MPI_STATUS_IGNORE);
  if ((my_rank/ncols)%2 == 1) {                         /*  send up-down                          */
    MPI_Irecv(p_a, 1, row, my_rank-ncols, 0, MPI_COMM_WORLD, &recv[4]);
    MPI_Isend(p_e, 1, row, my_rank-ncols, 0, MPI_COMM_WORLD, &send[4]);
    if (my_rank/ncols != nrows - 1) {
      MPI_Isend(p_f, 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, &send[5]);
      MPI_Irecv(p_g, 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, &recv[5]);
    }
  }
  else if (my_rank/ncols != nrows -1) {
    MPI_Isend(p_f, 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, &send[6]);
    MPI_Irecv(p_g, 1, row, my_rank+ncols, 0, MPI_COMM_WORLD, &recv[6]);
    if (my_rank/ncols != 0) {
      MPI_Irecv(p_a, 1, row, my_rank-ncols, 0, MPI_COMM_WORLD, &recv[7]);
      MPI_Isend(p_e, 1, row, my_rank-ncols, 0, MPI_COMM_WORLD, &send[7]);
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// updates the playing board
void update (int iteration){
  int x, y, xb, yb, neighbor;
  int *pointer_old = (iteration%2==0)?field_a:field_b;  /* pointer to the field that is current   */
  int *pointer_new = (iteration%2==0)?field_b:field_a;  /* pointer to the field that will update  */
  for (y = 0; y < local_height; y++) {                  /* loops through the local data           */
    for (x = 0; x < local_width; x++) {
      yb = y + 1;                                       /* shifts needed becuase the board        */
      xb = x + 1;                                       /*   is bigger due to ghost rows          */
      neighbor = 0;                                     /* initialize the neighbor count          */
      neighbor += pointer_old[(yb-1)*field_width+xb+1]; 
      neighbor += pointer_old[(yb-1)*field_width+xb]; 
      neighbor += pointer_old[(yb-1)*field_width+xb-1];
      neighbor += pointer_old[(yb)*field_width+xb+1];
      neighbor += pointer_old[(yb)*field_width+xb-1];
      neighbor += pointer_old[(yb+1)*field_width+xb+1];
      neighbor += pointer_old[(yb+1)*field_width+xb];
      neighbor += pointer_old[(yb+1)*field_width+xb-1];
      // sets the new filed to the old field value
      pointer_new[yb*field_width+xb]=pointer_old[yb*field_width+xb];
      // determines if the value of the new field should change
      if ((((neighbor < 2) || (neighbor > 3)) && (pointer_old[yb*field_width+xb])) || ((neighbor == 3) && (!pointer_old[yb*field_width+xb]))) {
        pointer_new[yb*field_width+xb] ^= 1;
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  char in_file[1000];
  char partition[100];
  int iterations, m_interval, w_interval, t_interval, offset, i;
  double start, finish, loc_elapsed, elapsed;
  MPI_Request send[8], recv[8];

  MPI_Init(&argc, &argv);                               /* start up MPI                           */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);               /* get the number of processes            */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes    */

  for (i = 0; i < 8; i++) {                             /* initializes the send and recv Request  */
    send[i] = MPI_REQUEST_NULL;                         /*  variables to NULL                     */
    recv[i] = MPI_REQUEST_NULL;
  }

  if (argc < 7) {                                       /* parses command line arguments          */
    cleanup(my_rank, "Error:  Too few arguments");
  }
  else if (argc == 7) {
    strcpy(in_file,   argv[1]);
    strcpy(partition, argv[2]);
    if (!((strcmp(partition, "slice") == 0) || (strcmp(partition, "checkerboard") == 0) || (strcmp(partition, "none") == 0))) {
      cleanup(my_rank, "Error:  Incorrect partition option.  Enter 'slice', 'checkerboard', or 'none'");
    }
    iterations = atoi(argv[3]);
    m_interval = atoi(argv[4]);
    w_interval = atoi(argv[5]);
    t_interval = atoi(argv[6]);
  }
  else if (argc > 7) {
    cleanup(my_rank, "Error:  Too many arguments");
  }

  fileread(in_file, partition, &offset);                /* function call to read the file in      */


  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  for (i = 0; i < iterations; i++) {
    summonspectre(i, send, recv);                       /* exchanges ghost fields                 */
    if (m_interval > 0) {                               /* checks if the bugs should be counted   */
      if (i%m_interval == 0) {
        measure(i);
      }
    }
    if ((w_interval > 0) && (i > 0)) {                  /* checks if the board should be written  */
      if (i%w_interval == 0) {
        filewrite(in_file, i, offset);
      }
    }

    if ((t_interval > 0) && (i > 0)) {                  /* measures the elapsed time              */
      if (i%t_interval == 0) {
        finish = MPI_Wtime();
        loc_elapsed = finish-start;
        MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (my_rank == 0) printf("Elapsed time = %e\n", elapsed);
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
      }
    }

    MPI_Waitall(8, recv, MPI_STATUS_IGNORE);
    update(i);                                          /* updates the game board                 */
    MPI_Waitall(8, send, MPI_STATUS_IGNORE);

  }
  cleanup(my_rank, "Thank you for playing!");           /* closes the program                     */
  return 0;
}

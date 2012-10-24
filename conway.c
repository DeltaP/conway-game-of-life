/*
 * Gregory Petropoulos
 *
 * This is my Conway's Game of Life program for the midterm
 *
 * To compile:  mpicc -g -Wall -std=c99 -o CGL conway.c -lm -llmpe -lmpe
 * To run:  mpiexec -n 1 ./CGL <input file name> <partition> <iteration> <m_interval> <w_interval>
 *          <> -> mandatory
 *          [] -> optional
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "mpe.h"
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
char header[100];
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Ends the program on an error and prints message to node 0
void cleanup (int my_rank, const char *message) {
  if (my_rank == 0) printf("%s\n",message);
  MPE_Stop_log();
  MPI_Finalize();                                       /* kills mpi                            */
  exit(0);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// reads in file
void fileread (char *filename, char *partition, int *offset) {
  if (strcmp(partition, "slice") == 0) {                /* determines the data decomposition    */
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

  MPI_File fh;                                          /* sets up MPI input                  */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  int newline_count = 0;
  char read;
  *offset = 0;

  do {                                                  /* reads in the header                */
    MPI_File_read_at_all(fh, *offset, &read, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    header[*offset]=read;
    if (read == '\n') newline_count++;
    (*offset)++;
    if (*offset == 100) cleanup(my_rank, "Error:  Header exceeds 100 characters, check file or recompile code");
  } while (newline_count < 3);
    
  char head[10];                                        /* parses the header and error      */
  int depth;                                            /*   checks                         */
  int rv               =  sscanf(header, "%6s\n%i %i\n%i\n", head, &width, &height, &depth);
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
  MPI_Offset disp = *offset;                            /* the initial displacement of      */
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

  int x, y, b, ll; 
  for (y = 0; y < field_height; y++ ) {             /* loops through the field          */
    for(x = 0; x < field_width; x++ ) {
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
// writes the game board out to file
void filewrite (char *in_file, int iteration, int offset) {
  char filename[1000];
  sprintf(filename, "%d_", iteration); 
  strcat(filename, in_file);

  MPI_File fh;                                          /* sets up MPI input                */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  if (my_rank == 0) {      /* writes the header                */
    MPI_File_write(fh,header,offset,MPI_CHAR,MPI_STATUS_IGNORE); 
  }

  char *temp=(char *)malloc( local_width * local_height * sizeof(char));
  int *field_pointer = (iteration%2==0) ? (field_a) : (field_b);

  int x, y, b, ll;
  for (y = 0; y < field_height; y++ ) {             /* loops through the field          */
    for(x = 0; x < field_width; x++ ) {
      if ((x != 0) && (y != 0) && (x != field_width - 1) && (y != field_height - 1)) {
        ll = (y * field_width + x);                     /* puts zeros at the boarders       */
        b = field_pointer[ll];
        b = (b==0)?0xFF:0;                              /* black = bugs; other = no bug     */
        temp[x-1+(y-1)*local_width] = (char)b;
      }
    }
  }

  MPI_Aint extent;                                      /* declares the extent              */
  MPI_Datatype etype, filetype, contig;                 /* derrived data types for IO       */
  MPI_Offset disp = offset;                             /* the initial displacement of      */

  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;

  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR       */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                  */
  extent = width * sizeof(char);                        /* total size of repeatable unit    */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data */
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


// -----------------------------------------------------------------
// swaps ghost rows --> looks messy is there a better way?
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
    MPI_Recv(pointer                              , 1, col, my_rank - 1     , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank%ncols != ncols -1) {             
    MPI_Send(local_width+pointer                  , 1, col, my_rank + 1     , 0, MPI_COMM_WORLD);
  }

  if ((my_rank%ncols)%2 == 1) {                         /* send left odd                    */
    MPI_Send(1+pointer                            , 1, col, my_rank-1       , 0, MPI_COMM_WORLD);
  }
  else if (my_rank%ncols != ncols -1) {
    MPI_Recv(local_width+1+pointer                , 1, col, my_rank+1       , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (my_rank%ncols != 0 && (my_rank%ncols)%2 == 0) {   /* send right even                  */
    MPI_Recv(pointer                              , 1, col, my_rank - 1     , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank%ncols != ncols -1 && (my_rank%ncols)%2==1) {
    MPI_Send(local_width+pointer                  , 1, col, my_rank + 1     , 0, MPI_COMM_WORLD);
  }

  if (my_rank%ncols != 0 && (my_rank%ncols)%2 == 0) {   /* send left even                   */
    MPI_Send(1+pointer                            , 1, col, my_rank-1       , 0, MPI_COMM_WORLD);
  }
  else if (my_rank%ncols != ncols -1 && (my_rank%ncols)%2==1) {
    MPI_Recv(local_width+1+pointer                , 1, col, my_rank+1       , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if ((my_rank/ncols)%2 == 1) {                         /* send down odd                    */
    MPI_Recv(pointer                              , 1, row, my_rank - ncols , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank/ncols != nrows -1) {
    MPI_Send(field_width*local_height+pointer     , 1, row, my_rank + ncols , 0, MPI_COMM_WORLD);
  }

  if ((my_rank/ncols)%2 == 1) {                         /* send up odd                      */
    MPI_Send(field_width+pointer                  , 1, row, my_rank-ncols   , 0, MPI_COMM_WORLD);
  }
  else if (my_rank/ncols != nrows -1) {
    MPI_Recv(field_width*(local_height+1)+pointer , 1, row, my_rank+ncols   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (my_rank/ncols != 0 && (my_rank/ncols)%2 == 0) {   /* send down even                   */
    MPI_Recv(pointer                              , 1, row, my_rank - ncols , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else if (my_rank/ncols != nrows -1 && (my_rank/ncols)%2==1) {
    MPI_Send(field_width*local_height+pointer     , 1, row, my_rank + ncols , 0, MPI_COMM_WORLD);
  }

  if (my_rank/ncols != 0 && (my_rank/ncols)%2 == 0) {   /* send up even                     */
    MPI_Send(field_width+pointer                  , 1, row, my_rank-ncols   , 0, MPI_COMM_WORLD);
  }
  else if (my_rank/ncols != nrows -1 && (my_rank/ncols)%2==1) {
    MPI_Recv(field_width*(local_height+1)+pointer , 1, row, my_rank+ncols   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// the main program
int main(int argc, char *argv[]) {
  char in_file[1000];
  char partition[100];
  int iterations, m_interval, w_interval, neighbor, offset, i, x, y;

  MPI_Init(&argc, &argv);                               /* start up MPI                         */
  MPE_Start_log();
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);               /* get the number of processes          */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes  */

  int event1a = MPE_Log_get_event_number();             /* creates events that are profiled     */
  int event1b = MPE_Log_get_event_number();
  int event2a = MPE_Log_get_event_number();
  int event2b = MPE_Log_get_event_number();
  int event3a = MPE_Log_get_event_number();
  int event3b = MPE_Log_get_event_number();
  int event4a = MPE_Log_get_event_number();
  int event4b = MPE_Log_get_event_number();
  int event5a = MPE_Log_get_event_number();
  int event5b = MPE_Log_get_event_number();

  if (my_rank == 0) {                                   /* gives descriptions to profiled events*/
    MPE_Describe_state(event1a, event1b, "Load File", "gray");
    MPE_Describe_state(event2a, event2b, "Exchange Ghost Rows", "green");
    MPE_Describe_state(event3a, event3b, "Measure", "red");
    MPE_Describe_state(event4a, event4b, "Write", "yellow");
    MPE_Describe_state(event5a, event5b, "Update", "blue");
  }

  if (argc < 6) {                                       /* parses command line arguments        */
    cleanup(my_rank, "Error:  Too few arguments");
  }
  else if (argc == 6) {
    strcpy(in_file,   argv[1]);
    strcpy(partition, argv[2]);
    if (!((strcmp(partition, "slice") == 0) || (strcmp(partition, "checkerboard") == 0) || (strcmp(partition, "none") == 0))) {
      cleanup(my_rank, "Error:  Incorrect partition option.  Enter 'slice', 'checkerboard', or 'none'");
    }
    iterations = atoi(argv[3]);
    m_interval = atoi(argv[4]);
    w_interval = atoi(argv[5]);
  }
  else if (argc > 6) {
    cleanup(my_rank, "Error:  Too many arguments");
  }

  MPE_Log_event(event1a, 0, "start Read");
  fileread(in_file, partition, &offset);                /* function call to read the file in    */
  MPE_Log_event(event1b, 0, "end Read");

  for (i = 0; i < iterations; i++) {

    MPE_Log_event(event2a, 0, "start Exchange Ghost Rows");
    summonspectre(i);                                   /* exchanges ghost fields               */
    MPE_Log_event(event2b, 0, "end Exchange Ghost Rows");

    if (m_interval > 0) {
      if (i%m_interval == 0) {
        MPE_Log_event(event3a, 0, "start Measure");
        measure(i);
        MPE_Log_event(event3b, 0, "end Measure");
      }
    }
    if (w_interval > 0) {
      if (i%w_interval == 0) {
        MPE_Log_event(event4a, 0, "start Write");
        filewrite(in_file, i, offset);
        MPE_Log_event(event4b, 0, "end Write");
      }
    }

    MPE_Log_event(event5a, 0, "start Update");
    int *pointer_old = (i%2==0)?field_a:field_b;        /* pointer to the field that is current */
    int *pointer_new = (i%2==0)?field_b:field_a;        /* pointer to the field that will update*/

    for (y = 0; y < local_height; y++) {            /* loops through the local data         */
      for (x = 0; x < local_width; x++) {
        int yb = y+1;                                   /* shifts needed becuase the            */
        int xb = x+1;                                   /*   is bigger due to ghost rows        */
        neighbor = 0;                                   /* initialize the neighbor count        */
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
    MPE_Log_event(event5b, 0, "end Update");
  }
  cleanup(my_rank, "Thank you for playing!");         /* closes the program                     */
  return 0;
}

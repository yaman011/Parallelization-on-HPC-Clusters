#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#include "fdtd.h"

#include <mpi.h>

#define DEFAULT_MPI_ERROR 1001

// MPI global variables
int root = 0;
int cart_rank;
int dims[3] = {0, 0, 0};
int coords[3];
int neighbors[6];
MPI_Comm cart_comm;

// Allocate and reuse communication buffers
static double *send_buffer_x = NULL, *recv_buffer_x = NULL;
static double *send_buffer_y = NULL, *recv_buffer_y = NULL;
static double *send_buffer_z = NULL, *recv_buffer_z = NULL;

static void initialize_mpi(int argc, char *argv[]) {
  // MPI local variables
  int reorder = 0;
  int world_size;
  int world_rank;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  int periods[3] = {0, 0, 0};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Get_processor_name(processor_name, &name_len);

  MPI_Dims_create(world_size, 3, dims);

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &cart_comm);
  MPI_Comm_rank(cart_comm, &cart_rank);

  MPI_Cart_coords(cart_comm, cart_rank, 3, coords);

  MPI_Cart_shift(cart_comm, 0, 1, &neighbors[LEFT], &neighbors[RIGHT]);
  MPI_Cart_shift(cart_comm, 1, 1, &neighbors[FRONT], &neighbors[BACK]);
  MPI_Cart_shift(cart_comm, 2, 1, &neighbors[DOWN], &neighbors[UP]);
}

void finalize_mpi() { MPI_Finalize(); }

void abort_mpi(int errorcode) { MPI_Abort(MPI_COMM_WORLD, errorcode); }

int main(int argc, char *argv[]) {
  // Initialize MPI
  initialize_mpi(argc, argv);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Print the number of OpenMP threads
  if (world_rank == 0) {
        printf("Number of OpenMP threads: %d\n", omp_get_max_threads());
    }
  if (argc < 2) {
    if (cart_rank == root) {
      printf("\nUsage: ./fdtd <param_file>\n\n");
    }
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  // Simulation and grid initialization
  global_grid_t global_grid;
  simulation_data_t simdata;
  init_simulation(&simdata, &global_grid, argv[1]);
  

  int gsize;
  MPI_Comm_size(cart_comm, &gsize);
  int rank_coords[3];
  int *displs = NULL, *rcounts = NULL;
  double *rbuff = NULL;
  global_data_t *global_output_data = NULL;

  int gnumnodesx = global_grid.numnodesx;
  int gnumnodesy = global_grid.numnodesy;
  int gnumnodesz = global_grid.numnodesz;
  int numnodestot = gnumnodesx * gnumnodesy * gnumnodesz;

  // Memory allocation for root
  if (cart_rank == root) {
    rbuff = (double *)malloc(numnodestot * sizeof(double));
    displs = (int *)malloc(gsize * sizeof(int));
    rcounts = (int *)malloc(gsize * sizeof(int));
    global_output_data = allocate_global_data(&global_grid);

    if (!rbuff || !displs || !rcounts || !global_output_data) {
      fprintf(stderr, "Memory allocation failed. Aborting...\n");
      abort_mpi(DEFAULT_MPI_ERROR);
    }

    for (int i = 0; i < gsize; i++) {
      MPI_Cart_coords(cart_comm, i, 3, rank_coords);

      // X
      int start_x = gnumnodesx * rank_coords[0] / dims[0];
      int end_x = gnumnodesx * (rank_coords[0] + 1) / dims[0] - 1;
      int lnumnodesx = (end_x - start_x + 1);
      rcounts[i] = lnumnodesx;

      // Y
      int start_y = gnumnodesy * rank_coords[1] / dims[1];
      int end_y = gnumnodesy * (rank_coords[1] + 1) / dims[1] - 1;
      int lnumnodesy = (end_y - start_y + 1);
      rcounts[i] *= lnumnodesy;

      // Z
      int start_z = gnumnodesz * rank_coords[2] / dims[2];
      int end_z = gnumnodesz * (rank_coords[2] + 1) / dims[2] - 1;
      int lnumnodesz = (end_z - start_z + 1);
      rcounts[i] *= lnumnodesz;

      displs[i] = (i == 0) ? 0 : displs[i - 1] + rcounts[i - 1];
    }
  }

  int numtimesteps = floor(simdata.params.maxt / simdata.params.dt);
  double start = GET_TIME();

  for (int tstep = 0; tstep <= numtimesteps; tstep++) {
    apply_source(&simdata, &global_grid, tstep);

    if (simdata.params.outrate > 0 && (tstep % simdata.params.outrate) == 0) {
      for (int i = 0; i < simdata.params.numoutputs; i++) {
        local_data_t *local_output_data = NULL;

        switch (simdata.params.outputs[i].source) {
        case PRESSURE:
          local_output_data = simdata.pold;
          break;
        case VELOCITYX:
          local_output_data = simdata.vxold;
          break;
        case VELOCITYY:
          local_output_data = simdata.vyold;
          break;
        case VELOCITYZ:
          local_output_data = simdata.vzold;
          break;

        default:
          break;
        }

        int lnumnodestot = NUMNODESTOT(simdata.pold->grid);
        MPI_Gatherv(local_output_data->vals, lnumnodestot, MPI_DOUBLE, rbuff,
                    rcounts, displs, MPI_DOUBLE, root, cart_comm);

        if (cart_rank == root) {
          for (int i = 0; i < gsize; i++) {
            MPI_Cart_coords(cart_comm, i, 3, rank_coords);

            // X
            int start_x = gnumnodesx * rank_coords[0] / dims[0];
            int end_x = gnumnodesx * (rank_coords[0] + 1) / dims[0] - 1;
            int lnumnodesx = (end_x - start_x + 1);

            // Y
            int start_y = gnumnodesy * rank_coords[1] / dims[1];
            int end_y = gnumnodesy * (rank_coords[1] + 1) / dims[1] - 1;
            int lnumnodesy = (end_y - start_y + 1);

            // Z
            int start_z = gnumnodesz * rank_coords[2] / dims[2];
            int end_z = gnumnodesz * (rank_coords[2] + 1) / dims[2] - 1;
            int lnumnodesz = (end_z - start_z + 1);

            for (int p = 0; p < lnumnodesz; p++) {
              for (int n = 0; n < lnumnodesy; n++) {
                for (int m = 0; m < lnumnodesx; m++) {
                  SETVALUE(global_output_data, start_x + m, start_y + n,
                           start_z + p,
                           rbuff[displs[i] + p * lnumnodesy * lnumnodesx +
                                 n * lnumnodesx + m]);
                }
              }
            }
          }

          double time = tstep * simdata.params.dt;
          write_output(&simdata.params.outputs[i], global_output_data, tstep,
                       time);
        }
      }
    }

    // Output progress
    if (cart_rank == root && tstep > 0 && tstep % (numtimesteps / 10) == 0) {
      printf("step %8d/%d", tstep, numtimesteps);

      if (tstep != numtimesteps) {
        double elapsed_sofar = GET_TIME() - start;
        double timeperstep_sofar = elapsed_sofar / tstep;

        double eta = (numtimesteps - tstep) * timeperstep_sofar;

        printf(" (ETA: %8.3lf seconds)", eta);
      }

      printf("\n");
      fflush(stdout);
    }

    // Update simulation state
    update_pressure(&simdata);
    update_velocities(&simdata, &global_grid);
    swap_timesteps(&simdata);
  }

  if (cart_rank == root) {
    double elapsed = GET_TIME() - start;
    double numupdates =
        (double)NUMNODESTOT(global_output_data->grid) * (numtimesteps + 1);
    double updatespers = numupdates / elapsed / 1e6;
    printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed,
           updatespers);
    fflush(stdout);

    free(rbuff);
    free(displs);
    free(rcounts);
    free(global_output_data->vals);
    free(global_output_data);
  }

  finalize_simulation(&simdata);
  finalize_mpi();

  return 0;
}

/******************************************************************************
 * Utilities functions                                                        *
 ******************************************************************************/

char *copy_string(char *str) {
  size_t len;
  if (str == NULL || (len = strlen(str)) == 0) {
    DEBUG_PRINT("NULL of zero length string passed as argument");
    return NULL;
  }

  char *cpy;
  if ((cpy = malloc((len + 1) * sizeof(char))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return NULL;
  }

  return strcpy(cpy, str);
}

void closest_index(global_grid_t *grid, double x, double y, double z, int *cx,
                   int *cy, int *cz) {
  int m = (int)((x - grid->xmin) / (grid->xmax - grid->xmin) * grid->numnodesx);
  int n = (int)((y - grid->ymin) / (grid->ymax - grid->ymin) * grid->numnodesy);
  int p = (int)((z - grid->zmin) / (grid->zmax - grid->zmin) * grid->numnodesz);

  *cx = (m < 0) ? 0 : (m > grid->numnodesx - 1) ? grid->numnodesx - 1 : m;
  *cy = (n < 0) ? 0 : (n > grid->numnodesy - 1) ? grid->numnodesy - 1 : n;
  *cz = (p < 0) ? 0 : (p > grid->numnodesz - 1) ? grid->numnodesz - 1 : p;
}

void lattice_below(global_grid_t *grid, double x, double y, double z, int *cx,
                   int *cy, int *cz) {
  int m = (int)(floor(x));
  int n = (int)(floor(y));
  int p = (int)(floor(z));

  *cx = (m < 0) ? 0 : (m > grid->numnodesx - 1) ? grid->numnodesx - 1 : m;
  *cy = (n < 0) ? 0 : (n > grid->numnodesy - 1) ? grid->numnodesy - 1 : n;
  *cz = (p < 0) ? 0 : (p > grid->numnodesz - 1) ? grid->numnodesz - 1 : p;
}

void print_source(source_t *source) {
  printf(" Source infos:\n\n");

  if (source->type == AUDIO) {
    double duration = (double)source->numsamples / source->sampling;

    printf("          type: audio data file\n");
    printf("      sampling: %d Hz\n", source->sampling);
    printf("      duration: %g\n", duration);

  } else {
    printf("          type: sine wave\n");
    printf("     frequency: %g Hz\n", source->data[0]);
  }

  printf("    position x: %g\n", source->posx);
  printf("    position y: %g\n", source->posy);
  printf("    position z: %g\n\n", source->posy);
}

void print_output(output_t *output) {
  switch (output->source) {
  case PRESSURE:
    printf("      pressure: ");
    break;
  case VELOCITYX:
    printf("    velocity X: ");
    break;
  case VELOCITYY:
    printf("    velocity Y: ");
    break;
  case VELOCITYZ:
    printf("    velocity Z: ");
    break;

  default:
    break;
  }

  switch (output->type) {
  case ALL:
    printf("complete dump");
    break;
  case CUTX:
    printf("cut along the x axis at %g", output->posx);
    break;
  case CUTY:
    printf("cut along the y axis at %g", output->posy);
    break;
  case CUTZ:
    printf("cut along the z axis at %g", output->posz);
    break;
  case POINT:
    printf("single point at %g %g %g", output->posx, output->posy,
           output->posz);
    break;

  default:
    break;
  }

  printf(" to file %s\n", output->filename);
}

/******************************************************************************
 * Data functions                                                             *
 ******************************************************************************/
global_data_t *allocate_global_data(global_grid_t *grid) {
  size_t numnodes = NUMNODESTOT(*grid);
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return NULL;
  }

  global_data_t *data;
  if ((data = malloc(sizeof(global_data_t))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data);
    return NULL;
  }

  if ((data->vals = malloc(numnodes * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data->vals);
    free(data);
    return NULL;
  }

  data->grid = *grid;

  return data;
}

local_data_t *allocate_local_data(local_grid_t *grid) {
  size_t numnodes = NUMNODESTOT(*grid);

  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return NULL;
  }

  local_data_t *data;
  if ((data = malloc(sizeof(local_data_t))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data);
    return NULL;
  }

  if ((data->vals = malloc(numnodes * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data->vals);
    free(data);
    return NULL;
  }

  data->grid = *grid;

  return data;
}

void fill_local_data(local_data_t *data, double value) {
  if (data == NULL) {
    DEBUG_PRINT("Invalid NULL data");
    return;
  }

  for (int m = 0; m < NUMNODESX(data); m++) {
    for (int n = 0; n < NUMNODESY(data); n++) {
      for (int p = 0; p < NUMNODESZ(data); p++) {
        SETVALUE(data, m, n, p, value);
      }
    }
  }
}

/******************************************************************************
 * Data file functions                                                        *
 ******************************************************************************/

FILE *create_datafile(global_grid_t grid, char *filename) {
  if (filename == NULL) {
    DEBUG_PRINT("Invalid NULL filename");
    return NULL;
  }

  FILE *fp;
  if ((fp = fopen(filename, "wb")) == NULL) {
    DEBUG_PRINTF("Failed to open file '%s'", filename);
    return NULL;
  }

  if (fwrite(&grid.numnodesx, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.numnodesy, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.numnodesz, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.xmin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.xmax, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.ymin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.ymax, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.zmin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.zmax, sizeof(double), 1, fp) != 1) {

    DEBUG_PRINTF("Failed to write header of file '%s'", filename);
    fclose(fp);
    return NULL;
  }

  return fp;
}

FILE *open_datafile(global_grid_t *grid, int *numsteps, char *filename) {
  if (grid == NULL || filename == NULL) {
    DEBUG_PRINT("Invalid NULL grid or filename");
    return NULL;
  }

  FILE *fp;
  if ((fp = fopen(filename, "rb")) == NULL) {
    DEBUG_PRINTF("Failed to open file '%s'", filename);
    return NULL;
  }

  fseek(fp, 0, SEEK_END);
  size_t file_size = ftell(fp);
  rewind(fp);

  if (fread(&grid->numnodesx, sizeof(int), 1, fp) != 1 ||
      fread(&grid->numnodesy, sizeof(int), 1, fp) != 1 ||
      fread(&grid->numnodesz, sizeof(int), 1, fp) != 1 ||
      fread(&grid->xmin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->xmax, sizeof(double), 1, fp) != 1 ||
      fread(&grid->ymin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->ymax, sizeof(double), 1, fp) != 1 ||
      fread(&grid->zmin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->zmax, sizeof(double), 1, fp) != 1) {
    DEBUG_PRINTF("Failed to read header of file '%s'", filename);
    fclose(fp);
    return NULL;
  }

  size_t numnodestot =
      (size_t)grid->numnodesx * grid->numnodesy * grid->numnodesz;

  size_t values_size = numnodestot * sizeof(double);
  size_t stepindex_size = sizeof(int);
  size_t timestamp_size = sizeof(double);
  size_t header_size = 6 * sizeof(double) + 3 * sizeof(int);

  size_t onetimestep_size = values_size + stepindex_size + timestamp_size;
  size_t alltimestep_size = file_size - header_size;

  if (alltimestep_size % onetimestep_size != 0) {
    DEBUG_PRINTF("Data size is inconsistent with number of nodes (%lu, %lu)",
                 alltimestep_size, onetimestep_size);

    fclose(fp);
    return NULL;
  }

  if (numsteps != NULL) {
    *numsteps = (alltimestep_size / onetimestep_size);
  }

  return fp;
}

global_data_t *read_data(FILE *fp, global_grid_t *grid, int *step,
                         double *time) {
  if (fp == NULL) {
    DEBUG_PRINT("Invalid NULL file pointer");
    return NULL;
  }

  double ltime;
  int lstep;

  size_t numnodes = NUMNODESTOT(*grid);

  global_data_t *data;
  if ((data = allocate_global_data(grid)) == NULL) {
    DEBUG_PRINT("Failed to allocate data");
    return NULL;
  }

  if (fread(&lstep, sizeof(int), 1, fp) != 1 ||
      fread(&ltime, sizeof(double), 1, fp) != 1 ||
      fread(data->vals, sizeof(double), numnodes, fp) != numnodes) {
    DEBUG_PRINT("Failed to read data");
    free(data);
    return NULL;
  }

  if (step != NULL)
    *step = lstep;
  if (time != NULL)
    *time = ltime;

  return data;
}

int write_data(FILE *fp, global_data_t *data, int step, double time) {
  if (fp == NULL || data == NULL || data->vals == NULL) {
    DEBUG_PRINT("Invalid NULL data or file pointer");
    return 1;
  }

  size_t numnodes = NUMNODESTOT(data->grid);
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return 1;
  }

  if (fwrite(&step, sizeof(int), 1, fp) != 1 ||
      fwrite(&time, sizeof(double), 1, fp) != 1 ||
      fwrite(data->vals, sizeof(double), numnodes, fp) != numnodes) {
    DEBUG_PRINT("Failed to write data");
    return 1;
  }

  return 0;
}

/******************************************************************************
 * Output file functions                                                      *
 ******************************************************************************/

int write_output(output_t *output, global_data_t *data, int step, double time) {
  if (output == NULL || data == NULL) {
    DEBUG_PRINT("NULL pointer passed as argument");
    return 1;
  }

  output_type_t type = output->type;

  if (type == ALL) {
    return write_data(output->fp, data, step, time);
  }

  int m, n, p;
  closest_index(&data->grid, output->posx, output->posy, output->posz, &m, &n,
                &p);

  int startm = (type == CUTX || type == POINT) ? m : 0;
  int startn = (type == CUTY || type == POINT) ? n : 0;
  int startp = (type == CUTZ || type == POINT) ? p : 0;

  int endm = (type == CUTX || type == POINT) ? m + 1 : NUMNODESX(data);
  int endn = (type == CUTY || type == POINT) ? n + 1 : NUMNODESY(data);
  int endp = (type == CUTZ || type == POINT) ? p + 1 : NUMNODESZ(data);

  global_data_t *tmpdata = allocate_global_data(&output->grid);

  for (m = startm; m < endm; m++) {
    for (n = startn; n < endn; n++) {
      for (p = startp; p < endp; p++) {
        int tmpm = m - startm;
        int tmpn = n - startn;
        int tmpp = p - startp;

        SETVALUE(tmpdata, tmpm, tmpn, tmpp, GETVALUE(data, m, n, p));
      }
    }
  }

  int writeok = (write_data(output->fp, tmpdata, step, time) == 0);

  free(tmpdata->vals);
  free(tmpdata);

  if (writeok == 0) {
    DEBUG_PRINT("Failed to write output data");
    return 1;
  }

  return 0;
}

int open_outputfile(output_t *output, global_grid_t *simgrid) {
  if (output == NULL || simgrid == NULL) {
    DEBUG_PRINT("Invalid NULL pointer in argment");
    return 1;
  }

  global_grid_t grid;

  output_type_t type = output->type;

  grid.numnodesx = (type == POINT || type == CUTX) ? 1 : simgrid->numnodesx;
  grid.numnodesy = (type == POINT || type == CUTY) ? 1 : simgrid->numnodesy;
  grid.numnodesz = (type == POINT || type == CUTZ) ? 1 : simgrid->numnodesz;

  grid.xmin = (type == POINT || type == CUTX) ? output->posx : simgrid->xmin;
  grid.xmax = (type == POINT || type == CUTX) ? output->posx : simgrid->xmax;

  grid.ymin = (type == POINT || type == CUTY) ? output->posy : simgrid->ymin;
  grid.ymax = (type == POINT || type == CUTY) ? output->posy : simgrid->ymax;

  grid.zmin = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmin;
  grid.zmax = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmax;

  FILE *fp;
  if ((fp = create_datafile(grid, output->filename)) == NULL) {
    DEBUG_PRINTF("Failed to open output file: '%s'", output->filename);
    return 1;
  }

  output->grid = grid;
  output->fp = fp;

  return 0;
}

/******************************************************************************
 * Parameter file functions                                                   *
 ******************************************************************************/

int read_audiosource(char *filename, source_t *source) {
  FILE *fp;
  if ((fp = fopen(filename, "rb")) == NULL) {
    DEBUG_PRINTF("Could not open source file '%s'", filename);
    return 1;
  }

  fseek(fp, 0, SEEK_END);
  size_t filesize = ftell(fp);
  rewind(fp);

  int numsamples = (filesize - sizeof(int)) / sizeof(double);

  int sampling;
  if (fread(&sampling, sizeof(int), 1, fp) != 1) {
    DEBUG_PRINT("Failed to read source data");
    fclose(fp);
    return 1;
  }

  double *data;
  if ((data = malloc(numsamples * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory for source data");
    return 1;
  }

  int readok = (fread(data, sizeof(double), numsamples, fp) == numsamples);

  fclose(fp);

  if (readok == 0) {
    DEBUG_PRINT("Failed to read source data");
    return 1;
  }

  source->data = data;
  source->numsamples = numsamples;
  source->sampling = sampling;

  return 0;
}

int read_outputparam(FILE *fp, output_t *output) {
  if (fp == NULL || output == NULL) {
    DEBUG_PRINT("NULL passed as argement");
    return 1;
  }

  char typekeyword[BUFSZ_SMALL];
  char sourcekeyword[BUFSZ_SMALL];
  char filename[BUFSZ_LARGE];

  double posxyz[3] = {0.0, 0.0, 0.0};

  if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1 ||
      fscanf(fp, BUFFMT_SMALL, sourcekeyword) != 1 ||
      fscanf(fp, BUFFMT_LARGE, filename) != 1) {

    DEBUG_PRINT("Failed to read an output parameter");
    return 1;
  }

  output_type_t type = CUTX;
  while (type < OUTPUT_TYPE_END &&
         strcmp(output_type_keywords[type], typekeyword) != 0) {
    type++;
  }

  if (type == OUTPUT_TYPE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
    return 1;
  }

  output_source_t source = PRESSURE;
  while (source < OUTPUT_SOURCE_END &&
         strcmp(output_source_keywords[source], sourcekeyword) != 0) {
    source++;
  }

  if (source == OUTPUT_SOURCE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", sourcekeyword);
    return 1;
  }

  int readok = 1;
  switch (type) {
  case CUTX:
    readok = (fscanf(fp, "%lf", &posxyz[0]) == 1);
    break;
  case CUTY:
    readok = (fscanf(fp, "%lf", &posxyz[1]) == 1);
    break;
  case CUTZ:
    readok = (fscanf(fp, "%lf", &posxyz[2]) == 1);
    break;
  case ALL:
    break;

  case POINT:
    readok =
        (fscanf(fp, "%lf %lf %lf", &posxyz[0], &posxyz[1], &posxyz[2]) == 3);
    break;

  default:
    break;
  }

  if (readok == 0) {
    DEBUG_PRINT("Failed to read an output parameter");
    return 1;
  }

  output->filename = copy_string(filename);
  output->type = type;
  output->source = source;
  output->posx = posxyz[0];
  output->posy = posxyz[1];
  output->posz = posxyz[2];

  return 0;
}

int read_sourceparam(FILE *fp, source_t *source) {
  char typekeyword[BUFSZ_SMALL];
  char filename[BUFSZ_LARGE];

  double freq, posx, posy, posz;

  if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1) {
    DEBUG_PRINT("Failed to read the source parameter");
    return 1;
  }

  source_type_t type = SINE;
  while (type < SOURCE_TYPE_END &&
         strcmp(source_type_keywords[type], typekeyword) != 0) {
    type++;
  }

  if (type == SOURCE_TYPE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
    return 1;
  }

  int readok = 1;
  switch (type) {
  case SINE:
    readok = (fscanf(fp, "%lf", &freq) == 1);
    break;
  case AUDIO:
    readok = (fscanf(fp, BUFFMT_LARGE, filename) == 1);
    break;

  default:
    break;
  }

  if (readok == 0 || fscanf(fp, "%lf %lf %lf", &posx, &posy, &posz) != 3) {
    DEBUG_PRINT("Failed to read the source parameter");
    return 1;
  }

  switch (type) {
  case AUDIO:
    read_audiosource(filename, source);
    break;
  case SINE: {
    if ((source->data = malloc(sizeof(double))) == NULL) {
      DEBUG_PRINT("Failed to allocate memory");
      return 1;
    }

    source->data[0] = freq;
    source->numsamples = 1;

    break;
  }

  default:
    break;
  }

  source->type = type;
  source->posx = posx;
  source->posy = posy;
  source->posz = posz;

  return 0;
}

int read_paramfile(parameters_t *params, const char *filename) {
  if (params == NULL || filename == NULL) {
    DEBUG_PRINT("Invalid print_out params or filename");
    return 1;
  }

  int outrate, numoutputs = 0;

  double dx, dt, maxt;

  char cin_filename[BUFSZ_LARGE];
  char rhoin_filename[BUFSZ_LARGE];

  source_t source;
  output_t *outputs = NULL;

  if ((outputs = malloc(sizeof(output_t) * MAX_OUTPUTS)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  FILE *fp;
  if ((fp = fopen(filename, "r")) == NULL) {
    DEBUG_PRINTF("Could not open parameter file '%s'", filename);
    return 1;
  }

  int readok =
      ((fscanf(fp, "%lf", &dx) == 1) && (fscanf(fp, "%lf", &dt) == 1) &&
       (fscanf(fp, "%lf", &maxt) == 1) && (fscanf(fp, "%d", &outrate) == 1) &&
       (fscanf(fp, BUFFMT_LARGE, cin_filename) == 1) &&
       (fscanf(fp, BUFFMT_LARGE, rhoin_filename) == 1));

  readok = (readok != 0 && read_sourceparam(fp, &source) == 0 &&
            fscanf(fp, " ") == 0);

  while (readok != 0 && numoutputs < MAX_OUTPUTS && feof(fp) == 0) {
    readok = (read_outputparam(fp, &outputs[numoutputs++]) == 0 &&
              fscanf(fp, " ") == 0);
  }

  fclose(fp);

  if (readok == 0) {
    DEBUG_PRINT("Failed to read parameter file");
    free(outputs);
    return 1;
  }

  if (numoutputs == 0) {
    free(outputs);
    outputs = NULL;

  } else if ((outputs = realloc(outputs, sizeof(output_t) * numoutputs)) ==
             NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  params->dx = dx;
  params->dt = dt;
  params->maxt = maxt;
  params->outrate = outrate;
  params->cin_filename = copy_string(cin_filename);
  params->rhoin_filename = copy_string(rhoin_filename);
  params->source = source;
  params->numoutputs = numoutputs;
  params->outputs = outputs;

  return 0;
}

/******************************************************************************
 * Simulation related functions                                               *
 ******************************************************************************/

int interpolate_inputmaps(simulation_data_t *simdata, local_grid_t *simgrid,
                          global_data_t *cin, global_data_t *rhoin) {
  if (simdata == NULL || cin == NULL) {
    DEBUG_PRINT("Invalid NULL simdata or cin");
    return 1;
  }

  if ((simdata->c = allocate_local_data(simgrid)) == NULL ||
      (simdata->rho = allocate_local_data(simgrid)) == NULL ||
      (simdata->rhohalf = allocate_local_data(simgrid)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  double dx = simdata->params.dx;
  double dxd2 = simdata->params.dx / 2;

  // use trilinear interpolation to evaluate the local values of c and rho

  for (int m = 0; m < simgrid->numnodesx; m++) {
    for (int n = 0; n < simgrid->numnodesy; n++) {
      for (int p = 0; p < simgrid->numnodesz; p++) {

        // point we want to evaluate
        double x = XMIN(cin) + (m + simgrid->start_x) * dx;
        double y = YMIN(cin) + (n + simgrid->start_y) * dx;
        double z = ZMIN(cin) + (p + simgrid->start_z) * dx;

        // ---- velocity ----

        // closest lattice point below (c000)
        int mc0, nc0, pc0;
        lattice_below(&cin->grid, x, y, z, &mc0, &nc0, &pc0);

        // closest lattice point above (c111)
        int mc1 = mc0 + 1;
        int nc1 = nc0 + 1;
        int pc1 = pc0 + 1;

        double xd = (x - mc0) / (mc1 - mc0);
        double yd = (y - nc0) / (nc1 - nc0);
        double zd = (z - pc0) / (pc1 - pc0);

        double c000 = GETVALUE(cin, mc0, nc0, pc0);
        double c001 = GETVALUE(cin, mc0, nc0, pc1);
        double c010 = GETVALUE(cin, mc0, nc1, pc0);
        double c011 = GETVALUE(cin, mc0, nc1, pc1);
        double c100 = GETVALUE(cin, mc1, nc0, pc0);
        double c101 = GETVALUE(cin, mc1, nc0, pc1);
        double c110 = GETVALUE(cin, mc1, nc1, pc0);
        double c111 = GETVALUE(cin, mc1, nc1, pc1);

        double c00 = c000 * (1 - xd) + c100 * xd;
        double c01 = c001 * (1 - xd) + c101 * xd;
        double c10 = c010 * (1 - xd) + c110 * xd;
        double c11 = c011 * (1 - xd) + c111 * xd;

        double c0 = c00 * (1 - yd) + c10 * yd;
        double c1 = c01 * (1 - yd) + c11 * yd;

        double c = c0 * (1 - zd) + c1 * zd;

        SETVALUE(simdata->c, m, n, p, c);

        // ---- density ----

        // closest lattice point below
        lattice_below(&rhoin->grid, x, y, z, &mc0, &nc0, &pc0);

        // closest lattice point above
        mc1 = mc0 + 1;
        nc1 = nc0 + 1;
        pc1 = pc0 + 1;

        xd = (x - mc0) / (mc1 - mc0);
        yd = (y - nc0) / (nc1 - nc0);
        zd = (z - pc0) / (pc1 - pc0);

        double rho000 = GETVALUE(rhoin, mc0, nc0, pc0);
        double rho001 = GETVALUE(rhoin, mc0, nc0, pc1);
        double rho010 = GETVALUE(rhoin, mc0, nc1, pc0);
        double rho011 = GETVALUE(rhoin, mc0, nc1, pc1);
        double rho100 = GETVALUE(rhoin, mc1, nc0, pc0);
        double rho101 = GETVALUE(rhoin, mc1, nc0, pc1);
        double rho110 = GETVALUE(rhoin, mc1, nc1, pc0);
        double rho111 = GETVALUE(rhoin, mc1, nc1, pc1);

        double rho00 = rho000 * (1 - xd) + rho100 * xd;
        double rho01 = rho001 * (1 - xd) + rho101 * xd;
        double rho10 = rho010 * (1 - xd) + rho110 * xd;
        double rho11 = rho011 * (1 - xd) + rho111 * xd;

        double rho0 = rho00 * (1 - yd) + rho10 * yd;
        double rho1 = rho01 * (1 - yd) + rho11 * yd;

        double rho = rho0 * (1 - zd) + rho1 * zd;

        SETVALUE(simdata->rho, m, n, p, rho);

        // ---- density at half step ----

        x += dxd2;
        y += dxd2;
        z += dxd2;

        // closest lattice point below
        lattice_below(&rhoin->grid, x, y, z, &mc0, &nc0, &pc0);

        // closest lattice point above
        mc1 = mc0 + 1;
        nc1 = nc0 + 1;
        pc1 = pc0 + 1;

        xd = (x - mc0) / (mc1 - mc0);
        yd = (y - nc0) / (nc1 - nc0);
        zd = (z - pc0) / (pc1 - pc0);

        rho000 = GETVALUE(rhoin, mc0, nc0, pc0);
        rho001 = GETVALUE(rhoin, mc0, nc0, pc1);
        rho010 = GETVALUE(rhoin, mc0, nc1, pc0);
        rho011 = GETVALUE(rhoin, mc0, nc1, pc1);
        rho100 = GETVALUE(rhoin, mc1, nc0, pc0);
        rho101 = GETVALUE(rhoin, mc1, nc0, pc1);
        rho110 = GETVALUE(rhoin, mc1, nc1, pc0);
        rho111 = GETVALUE(rhoin, mc1, nc1, pc1);

        rho00 = rho000 * (1 - xd) + rho100 * xd;
        rho01 = rho001 * (1 - xd) + rho101 * xd;
        rho10 = rho010 * (1 - xd) + rho110 * xd;
        rho11 = rho011 * (1 - xd) + rho111 * xd;

        rho0 = rho00 * (1 - yd) + rho10 * yd;
        rho1 = rho01 * (1 - yd) + rho11 * yd;

        rho = rho0 * (1 - zd) + rho1 * zd;

        SETVALUE(simdata->rhohalf, m, n, p, rho);
      }
    }
  }

  return 0;
}

void apply_source(simulation_data_t *simdata, global_grid_t *global_grid,
                  int step) {
  source_t *source = &simdata->params.source;

  double posx = source->posx;
  double posy = source->posy;
  double posz = source->posz;

  double t = step * simdata->params.dt;

  int m, n, p;
  closest_index(global_grid, posx, posy, posz, &m, &n, &p);

  m = m - simdata->pold->grid.start_x;
  n = n - simdata->pold->grid.start_y;
  p = p - simdata->pold->grid.start_z;

  if (m < 0 || m > simdata->pold->grid.numnodesx - 1 || n < 0 ||
      n > simdata->pold->grid.numnodesy - 1 || p < 0 ||
      p > simdata->pold->grid.numnodesz - 1) {
    return;
  }

  if (source->type == SINE) {
    double freq = source->data[0];

    SETVALUE(simdata->pold, m, n, p, sin(2 * M_PI * freq * t));

  } else if (source->type == AUDIO) {
    int sample = MIN((int)(t * source->sampling), source->numsamples - 1);

    SETVALUE(simdata->pold, m, n, p, simdata->params.source.data[sample]);
  }
}

// PRESSURE UPDATE

static void update_inner_pressure(simulation_data_t *simdata) {
    const double dtdx = simdata->params.dt / simdata->params.dx;
    const int numnodesx = NUMNODESX(simdata->pold);
    const int numnodesy = NUMNODESY(simdata->pold);
    const int numnodesz = NUMNODESZ(simdata->pold);
    
    #pragma omp parallel for collapse(3)
    for (int p = 1; p < numnodesz; p++) {
        for (int n = 1; n < numnodesy; n++) {
            for (int m = 1; m < numnodesx; m++) {
                const double c = GETVALUE(simdata->c, m, n, p);
                const double rho = GETVALUE(simdata->rho, m, n, p);
                const double rhoc2dtdx = rho * c * c * dtdx;

                const double dvx = GETVALUE(simdata->vxold, m, n, p) - GETVALUE(simdata->vxold, m - 1, n, p);
                const double dvy = GETVALUE(simdata->vyold, m, n, p) - GETVALUE(simdata->vyold, m, n - 1, p);
                const double dvz = GETVALUE(simdata->vzold, m, n, p) - GETVALUE(simdata->vzold, m, n, p - 1);

                const double divv = dvx + dvy + dvz;
                const double prev_p = GETVALUE(simdata->pold, m, n, p);

                const double pnew = prev_p - rhoc2dtdx * divv;

                SETVALUE(simdata->pnew, m, n, p, pnew);
            }
        }
    }
}


static void update_outer_pressure(simulation_data_t *simdata, double *recv_buffer_x, double *recv_buffer_y, double *recv_buffer_z) {
    const int numnodesx = NUMNODESX(simdata->pold);
    const int numnodesy = NUMNODESY(simdata->pold);
    const int numnodesz = NUMNODESZ(simdata->pold);
    const double dtdx = simdata->params.dt / simdata->params.dx;

  
  
      #pragma omp parallel
    {
        // Boundary at x = 0
        #pragma omp for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int n = 0; n < numnodesy; n++) {
            int m = 0;
            const double c = GETVALUE(simdata->c, m, n, p);
            const double rho = GETVALUE(simdata->rho, m, n, p);
            const double rhoc2dtdx = rho * c * c * dtdx;

            double dvx = GETVALUE(simdata->vxold, m, n, p);
            if (neighbors[LEFT] != MPI_PROC_NULL) {
                dvx -= recv_buffer_x[p * numnodesy + n];
            }

            double dvy = GETVALUE(simdata->vyold, m, n, p) - GETVALUE(simdata->vyold, m, n - 1, p);
            double dvz = GETVALUE(simdata->vzold, m, n, p) - GETVALUE(simdata->vzold, m, n, p - 1);

            const double prev_p = GETVALUE(simdata->pold, m, n, p);
            const double pnew = prev_p - rhoc2dtdx * (dvx + dvy + dvz);
            SETVALUE(simdata->pnew, m, n, p, pnew);
        }
    }

           // Boundary at y = 0
        #pragma omp for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int m = 0; m < numnodesx; m++) {
            int n = 0;
            const double c = GETVALUE(simdata->c, m, n, p);
            const double rho = GETVALUE(simdata->rho, m, n, p);
            const double rhoc2dtdx = rho * c * c * dtdx;

            double dvy = GETVALUE(simdata->vyold, m, n, p);
            if (neighbors[FRONT] != MPI_PROC_NULL) {
                dvy -= recv_buffer_y[p * numnodesx + m];
            }

            double dvx = GETVALUE(simdata->vxold, m, n, p) - GETVALUE(simdata->vxold, m - 1, n, p);
            double dvz = GETVALUE(simdata->vzold, m, n, p) - GETVALUE(simdata->vzold, m, n, p - 1);

            const double prev_p = GETVALUE(simdata->pold, m, n, p);
            const double pnew = prev_p - rhoc2dtdx * (dvx + dvy + dvz);
            SETVALUE(simdata->pnew, m, n, p, pnew);
        }
    }

           // Boundary at z = 0
        #pragma omp for collapse(2)
    for (int n = 0; n < numnodesy; n++) {
        for (int m = 0; m < numnodesx; m++) {
            int p = 0;
            const double c = GETVALUE(simdata->c, m, n, p);
            const double rho = GETVALUE(simdata->rho, m, n, p);
            const double rhoc2dtdx = rho * c * c * dtdx;

            double dvz = GETVALUE(simdata->vzold, m, n, p);
            if (neighbors[DOWN] != MPI_PROC_NULL) {
                dvz -= recv_buffer_z[n * numnodesx + m];
            }

            double dvx = GETVALUE(simdata->vxold, m, n, p) - GETVALUE(simdata->vxold, m - 1, n, p);
            double dvy = GETVALUE(simdata->vyold, m, n, p) - GETVALUE(simdata->vyold, m, n - 1, p);

            const double prev_p = GETVALUE(simdata->pold, m, n, p);
            const double pnew = prev_p - rhoc2dtdx * (dvx + dvy + dvz);
            SETVALUE(simdata->pnew, m, n, p, pnew);
        }
    }
}
}


void update_pressure(simulation_data_t *simdata) {
    const int numnodesx = NUMNODESX(simdata->pold);
    const int numnodesy = NUMNODESY(simdata->pold);
    const int numnodesz = NUMNODESZ(simdata->pold);

    // Allocate and initialize communication buffers
    allocate_comm_buffers(numnodesx, numnodesy, numnodesz);

    // Pack the necessary values into send buffers
    #pragma omp parallel for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int n = 0; n < numnodesy; n++) {
            send_buffer_x[p * numnodesy + n] = GETVALUE(simdata->vxold, numnodesx - 1, n, p);
        }
    }

     #pragma omp parallel for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int m = 0; m < numnodesx; m++) {
            send_buffer_y[p * numnodesx + m] = GETVALUE(simdata->vyold, m, numnodesy - 1, p);
        }
    }



    #pragma omp parallel for collapse(2)
    for (int n = 0; n < numnodesy; n++) {
        for (int m = 0; m < numnodesx; m++) {
            send_buffer_z[n * numnodesx + m] = GETVALUE(simdata->vzold, m, n, numnodesz - 1);
        }
    }

    // Posting non-blocking sends and receives
    MPI_Request request_recv[3], request_send[3];
    MPI_Irecv(recv_buffer_x, numnodesy * numnodesz, MPI_DOUBLE, neighbors[LEFT], 0, cart_comm, &request_recv[0]);
    MPI_Irecv(recv_buffer_y, numnodesx * numnodesz, MPI_DOUBLE, neighbors[FRONT], 1, cart_comm, &request_recv[1]);
    MPI_Irecv(recv_buffer_z, numnodesx * numnodesy, MPI_DOUBLE, neighbors[DOWN], 2, cart_comm, &request_recv[2]);

    MPI_Isend(send_buffer_x, numnodesy * numnodesz, MPI_DOUBLE, neighbors[RIGHT], 0, cart_comm, &request_send[0]);
    MPI_Isend(send_buffer_y, numnodesx * numnodesz, MPI_DOUBLE, neighbors[BACK], 1, cart_comm, &request_send[1]);
    MPI_Isend(send_buffer_z, numnodesx * numnodesy, MPI_DOUBLE, neighbors[UP], 2, cart_comm, &request_send[2]);

    // Update inner pressure field
    update_inner_pressure(simdata);

    // Wait for all communications to complete
    MPI_Waitall(3, request_recv, MPI_STATUSES_IGNORE);

    // Update outer pressure field
    update_outer_pressure(simdata, recv_buffer_x, recv_buffer_y, recv_buffer_z);

    MPI_Waitall(3, request_send, MPI_STATUSES_IGNORE);

}


// VELOCITY UPDATE

static void update_inner_velocities(simulation_data_t *simdata) {
    const double dtdx = simdata->params.dt / simdata->params.dx;
    const int numnodesx = NUMNODESX(simdata->vxold);
    const int numnodesy = NUMNODESY(simdata->vxold);
    const int numnodesz = NUMNODESZ(simdata->vxold);
    
    #pragma omp parallel for collapse(3)
    for (int p = 0; p < numnodesz - 1; p++) {
        for (int n = 0; n < numnodesy - 1; n++) {
            for (int m = 0; m < numnodesx - 1; m++) {
                const double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);
                const double p_mnq = GETVALUE(simdata->pnew, m, n, p);

                const double dpx = GETVALUE(simdata->pnew, m + 1, n, p) - p_mnq;
                const double dpy = GETVALUE(simdata->pnew, m, n + 1, p) - p_mnq;
                const double dpz = GETVALUE(simdata->pnew, m, n, p + 1) - p_mnq;

                const double vxnew = GETVALUE(simdata->vxold, m, n, p) - dtdxrho * dpx;
                const double vynew = GETVALUE(simdata->vyold, m, n, p) - dtdxrho * dpy;
                const double vznew = GETVALUE(simdata->vzold, m, n, p) - dtdxrho * dpz;

                SETVALUE(simdata->vxnew, m, n, p, vxnew);
                SETVALUE(simdata->vynew, m, n, p, vynew);
                SETVALUE(simdata->vznew, m, n, p, vznew);
            }
        }
    }
}


static void update_outer_velocities(simulation_data_t *simdata,
                                    global_grid_t *global_grid,
                                    double *recv_buffer_x,
                                    double *recv_buffer_y,
                                    double *recv_buffer_z) {
  const int numnodesx = NUMNODESX(simdata->pold);
  const int numnodesy = NUMNODESY(simdata->pold);
  const int numnodesz = NUMNODESZ(simdata->pold);
  const double dtdx = simdata->params.dt / simdata->params.dx;
  double dpx, dpy, dpz;

    #pragma omp parallel
    {
        #pragma omp for collapse(2)
  for (int p = 0; p < numnodesz; p++) {
    for (int n = 0; n < numnodesy; n++) {
      int m = numnodesx - 1; // Boundary at x = Lx
      const double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);
      const int numnodesx = NUMNODESX(simdata->pold);
      const int numnodesy = NUMNODESY(simdata->pold);
      const int numnodesz = NUMNODESZ(simdata->pold);

      double p_mnp = GETVALUE(simdata->pnew, m, n, p);

      int mp1 = MIN(global_grid->numnodesx - 1, m + 1);
      int np1 = MIN(global_grid->numnodesy - 1, n + 1);
      int pp1 = MIN(global_grid->numnodesz - 1, p + 1);

      if (mp1 < numnodesx) {
        dpx = GETVALUE(simdata->pnew, mp1, n, p) - p_mnp;
      } else {
        dpx = recv_buffer_x[p * numnodesy + n] - p_mnp;
      }

      if (np1 < numnodesy) {
        dpy = GETVALUE(simdata->pnew, m, np1, p) - p_mnp;
      } else {
        dpy = recv_buffer_y[p * numnodesx + m] - p_mnp;
      }

      if (pp1 < numnodesz) {
        dpz = GETVALUE(simdata->pnew, m, n, pp1) - p_mnp;
      } else {
        dpz = recv_buffer_z[n * numnodesx + m] - p_mnp;
      }

      double prev_vx = GETVALUE(simdata->vxold, m, n, p);
      double prev_vy = GETVALUE(simdata->vyold, m, n, p);
      double prev_vz = GETVALUE(simdata->vzold, m, n, p);

      SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
      SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
      SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
  

  #pragma omp for collapse(2)
  for (int p = 0; p < numnodesz; p++) {
    for (int m = 0; m < numnodesx; m++) {
      int n = numnodesy - 1; // Boundary at y = Ly
      const double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);
      const int numnodesx = NUMNODESX(simdata->pold);
      const int numnodesy = NUMNODESY(simdata->pold);
      const int numnodesz = NUMNODESZ(simdata->pold);

      double p_mnp = GETVALUE(simdata->pnew, m, n, p);

      int mp1 = MIN(global_grid->numnodesx - 1, m + 1);
      int np1 = MIN(global_grid->numnodesy - 1, n + 1);
      int pp1 = MIN(global_grid->numnodesz - 1, p + 1);

      if (mp1 < numnodesx) {
        dpx = GETVALUE(simdata->pnew, mp1, n, p) - p_mnp;
      } else {
        dpx = recv_buffer_x[p * numnodesy + n] - p_mnp;
      }

      if (np1 < numnodesy) {
        dpy = GETVALUE(simdata->pnew, m, np1, p) - p_mnp;
      } else {
        dpy = recv_buffer_y[p * numnodesx + m] - p_mnp;
      }

      if (pp1 < numnodesz) {
        dpz = GETVALUE(simdata->pnew, m, n, pp1) - p_mnp;
      } else {
        dpz = recv_buffer_z[n * numnodesx + m] - p_mnp;
      }

      double prev_vx = GETVALUE(simdata->vxold, m, n, p);
      double prev_vy = GETVALUE(simdata->vyold, m, n, p);
      double prev_vz = GETVALUE(simdata->vzold, m, n, p);

      SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
      SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
      SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
  


  #pragma omp for collapse(2)
  for (int n = 0; n < numnodesy; n++) {
    for (int m = 0; m < numnodesx; m++) {
      int p = numnodesz - 1; // Boundary at z = Lz
      const double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);
      const int numnodesx = NUMNODESX(simdata->pold);
      const int numnodesy = NUMNODESY(simdata->pold);
      const int numnodesz = NUMNODESZ(simdata->pold);

      double p_mnp = GETVALUE(simdata->pnew, m, n, p);
      
      int mp1 = MIN(global_grid->numnodesx - 1, m + 1);
      int np1 = MIN(global_grid->numnodesy - 1, n + 1);
      int pp1 = MIN(global_grid->numnodesz - 1, p + 1);

      if (mp1 < numnodesx) {
        dpx = GETVALUE(simdata->pnew, mp1, n, p) - p_mnp;
      } else {
        dpx = recv_buffer_x[p * numnodesy + n] - p_mnp;
      }

      if (np1 < numnodesy) {
        dpy = GETVALUE(simdata->pnew, m, np1, p) - p_mnp;
      } else {
        dpy = recv_buffer_y[p * numnodesx + m] - p_mnp;
      }

      if (pp1 < numnodesz) {
        dpz = GETVALUE(simdata->pnew, m, n, pp1) - p_mnp;
      } else {
        dpz = recv_buffer_z[n * numnodesx + m] - p_mnp;
      }

      double prev_vx = GETVALUE(simdata->vxold, m, n, p);
      double prev_vy = GETVALUE(simdata->vyold, m, n, p);
      double prev_vz = GETVALUE(simdata->vzold, m, n, p);

      SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
      SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
      SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
    }
  }
}
}

void update_velocities(simulation_data_t *simdata, global_grid_t *global_grid) {
    const int numnodesx = NUMNODESX(simdata->pold);
    const int numnodesy = NUMNODESY(simdata->pold);
    const int numnodesz = NUMNODESZ(simdata->pold);
    // Allocate and initialize communication buffers
    allocate_comm_buffers(numnodesx, numnodesy, numnodesz);

    // Pack the necessary values into send buffers
    #pragma omp parallel for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int n = 0; n < numnodesy; n++) {
            send_buffer_x[p * numnodesy + n] = GETVALUE(simdata->pnew, 0, n, p);
        }
    }
   
    #pragma omp parallel for collapse(2)
    for (int p = 0; p < numnodesz; p++) {
        for (int m = 0; m < numnodesx; m++) {
            send_buffer_y[p * numnodesx + m] = GETVALUE(simdata->pnew, m, 0, p);
        }
    }
    

    #pragma omp parallel for collapse(2)
    for (int n = 0; n < numnodesy; n++) {
        for (int m = 0; m < numnodesx; m++) {
            send_buffer_z[n * numnodesx + m] = GETVALUE(simdata->pnew, m, n, 0);
        }
    }

    // Posting non-blocking sends and receives
    MPI_Request request_recv[3], request_send[3];
    MPI_Irecv(recv_buffer_x, numnodesy * numnodesz, MPI_DOUBLE, neighbors[RIGHT], 0, cart_comm, &request_recv[0]);
    MPI_Irecv(recv_buffer_y, numnodesx * numnodesz, MPI_DOUBLE, neighbors[BACK], 1, cart_comm, &request_recv[1]);
    MPI_Irecv(recv_buffer_z, numnodesx * numnodesy, MPI_DOUBLE, neighbors[UP], 2, cart_comm, &request_recv[2]);

    MPI_Isend(send_buffer_x, numnodesy * numnodesz, MPI_DOUBLE, neighbors[LEFT], 0, cart_comm, &request_send[0]);
    MPI_Isend(send_buffer_y, numnodesx * numnodesz, MPI_DOUBLE, neighbors[FRONT], 1, cart_comm, &request_send[1]);
    MPI_Isend(send_buffer_z, numnodesx * numnodesy, MPI_DOUBLE, neighbors[DOWN], 2, cart_comm, &request_send[2]);

    // Update inner velocities
    update_inner_velocities(simdata);

    // Wait for all communications to complete
    MPI_Waitall(3, request_recv, MPI_STATUSES_IGNORE);

    // Update outer velocities
    update_outer_velocities(simdata, global_grid, recv_buffer_x, recv_buffer_y, recv_buffer_z);

    MPI_Waitall(3, request_send, MPI_STATUSES_IGNORE);
}

void allocate_comm_buffers(int numnodesx, int numnodesy, int numnodesz) {
    if (send_buffer_x == NULL) {
        send_buffer_x = malloc(numnodesy * numnodesz * sizeof(double));
        recv_buffer_x = malloc(numnodesy * numnodesz * sizeof(double));
        send_buffer_y = malloc(numnodesx * numnodesz * sizeof(double));
        recv_buffer_y = malloc(numnodesx * numnodesz * sizeof(double));
        send_buffer_z = malloc(numnodesx * numnodesy * sizeof(double));
        recv_buffer_z = malloc(numnodesx * numnodesy * sizeof(double));

        if (!send_buffer_x || !recv_buffer_x || !send_buffer_y || !recv_buffer_y || !send_buffer_z || !recv_buffer_z) {
            fprintf(stderr, "Failed to allocate communication buffers\n");
            abort_mpi(DEFAULT_MPI_ERROR);
        }
    }
}

void free_comm_buffers() {
    if (send_buffer_x != NULL) {
        free(send_buffer_x);
        free(recv_buffer_x);
        free(send_buffer_y);
        free(recv_buffer_y);
        free(send_buffer_z);
        free(recv_buffer_z);
        
        send_buffer_x = recv_buffer_x = NULL;
        send_buffer_y = recv_buffer_y = NULL;
        send_buffer_z = recv_buffer_z = NULL;
    }
}



void init_simulation(simulation_data_t *simdata, global_grid_t *global_grid,
                     const char *params_filename) {
  if (read_paramfile(&simdata->params, params_filename) != 0) {
    printf("Failed to read parameters. Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  global_grid_t rhoin_grid;
  global_grid_t cin_grid;
  local_grid_t sim_grid;

  int rho_numstep;
  int c_numstep;

  // Open data files for rho and c maps
  FILE *rhofp =
      open_datafile(&rhoin_grid, &rho_numstep, simdata->params.rhoin_filename);
  FILE *cfp =
      open_datafile(&cin_grid, &c_numstep, simdata->params.cin_filename);

  if (rhofp == NULL || rho_numstep <= 0) {
    printf("Failed to open the density map file. Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  if (cfp == NULL || c_numstep <= 0) {
    printf("Failed to open the speed map file. Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  // Check if the grids are compatible
  if (rhoin_grid.xmin != cin_grid.xmin || rhoin_grid.ymin != cin_grid.ymin ||
      rhoin_grid.zmin != cin_grid.zmin || rhoin_grid.xmax != cin_grid.xmax ||
      rhoin_grid.ymax != cin_grid.ymax || rhoin_grid.zmax != cin_grid.zmax) {
    printf("Grids for the density and speed are not the same. "
           "Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  // Read data from input files
  global_data_t *rho_map = read_data(rhofp, &rhoin_grid, NULL, NULL);
  global_data_t *c_map = read_data(cfp, &cin_grid, NULL, NULL);

  if (rho_map == NULL || c_map == NULL) {
    printf("Failed to read data from input maps. Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  fclose(rhofp);
  fclose(cfp);

  // Global grid
  double gxmin = rhoin_grid.xmin;
  double gxmax = rhoin_grid.xmax;
  double gymin = rhoin_grid.ymin;
  double gymax = rhoin_grid.ymax;
  double gzmin = rhoin_grid.zmin;
  double gzmax = rhoin_grid.zmax;

  int gnumnodesx = MAX(floor((gxmax - gxmin) / simdata->params.dx), 1);
  int gnumnodesy = MAX(floor((gymax - gymin) / simdata->params.dx), 1);
  int gnumnodesz = MAX(floor((gzmax - gzmin) / simdata->params.dx), 1);

  global_grid->xmin = gxmin;
  global_grid->xmax = gxmax;
  global_grid->ymin = gymin;
  global_grid->ymax = gymax;
  global_grid->zmin = gzmin;
  global_grid->zmax = gzmax;

  global_grid->numnodesx = gnumnodesx;
  global_grid->numnodesy = gnumnodesy;
  global_grid->numnodesz = gnumnodesz;

  // Local grid
  // X
  int start_x = gnumnodesx * coords[0] / dims[0];
  int end_x = gnumnodesx * (coords[0] + 1) / dims[0] - 1;
  int lnumnodesx = end_x - start_x + 1;

  sim_grid.start_x = start_x;
  sim_grid.end_x = end_x;
  sim_grid.numnodesx = lnumnodesx;

  // Y
  int start_y = gnumnodesy * coords[1] / dims[1];
  int end_y = gnumnodesy * (coords[1] + 1) / dims[1] - 1;
  int lnumnodesy = end_y - start_y + 1;

  sim_grid.start_y = start_y;
  sim_grid.end_y = end_y;
  sim_grid.numnodesy = lnumnodesy;

  // Z
  int start_z = gnumnodesz * coords[2] / dims[2];
  int end_z = gnumnodesz * (coords[2] + 1) / dims[2] - 1;
  int lnumnodesz = end_z - start_z + 1;

  sim_grid.start_z = start_z;
  sim_grid.end_z = end_z;
  sim_grid.numnodesz = lnumnodesz;

  if (interpolate_inputmaps(simdata, &sim_grid, c_map, rho_map) != 0) {
    printf("Error while converting input map to simulation grid. "
           "Aborting...\n\n");
    abort_mpi(DEFAULT_MPI_ERROR);
  }

  if (cart_rank == root) {
    if (simdata->params.outrate > 0 && simdata->params.outputs != NULL) {
      for (int i = 0; i < simdata->params.numoutputs; i++) {
        char *outfilei = simdata->params.outputs[i].filename;

        for (int j = 0; j < i; j++) {
          char *outfilej = simdata->params.outputs[j].filename;

          if (strcmp(outfilei, outfilej) == 0) {
            printf("Duplicate output file: '%s'. Aborting...\n\n", outfilei);
            abort_mpi(DEFAULT_MPI_ERROR);
          }
        }
      }

      for (int i = 0; i < simdata->params.numoutputs; i++) {
        output_t *output = &simdata->params.outputs[i];

        if (open_outputfile(output, global_grid) != 0) {
          printf("Failed to open output file: '%s'. Aborting...\n\n",
                 output->filename);
          abort_mpi(DEFAULT_MPI_ERROR);
          exit(1);
        }
      }
    }
  }

  if ((simdata->pold = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->pnew = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vxold = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vxnew = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vyold = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vynew = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vzold = allocate_local_data(&sim_grid)) == NULL ||
      (simdata->vznew = allocate_local_data(&sim_grid)) == NULL) {
    printf("Failed to allocate memory. Aborting...\n\n");

    free(simdata->pold);
    free(simdata->pnew);
    free(simdata->vxold);
    free(simdata->vxnew);
    free(simdata->vyold);
    free(simdata->vynew);
    free(simdata->vzold);
    free(simdata->vznew);

    abort_mpi(DEFAULT_MPI_ERROR);
  }

  fill_local_data(simdata->pold, 0.0);
  fill_local_data(simdata->pnew, 0.0);

  fill_local_data(simdata->vynew, 0.0);
  fill_local_data(simdata->vxold, 0.0);
  fill_local_data(simdata->vynew, 0.0);
  fill_local_data(simdata->vyold, 0.0);
  fill_local_data(simdata->vznew, 0.0);
  fill_local_data(simdata->vzold, 0.0);

  if (cart_rank == root) {
    printf("\n");
    printf(" Grid spacing: %g\n", simdata->params.dx);
    printf("  Grid size X: %d\n", sim_grid.numnodesx);
    printf("  Grid size Y: %d\n", sim_grid.numnodesy);
    printf("  Grid size Z: %d\n", sim_grid.numnodesz);
    printf("    Time step: %g\n", simdata->params.dt);
    printf(" Maximum time: %g\n\n", simdata->params.maxt);

    if (simdata->params.outrate > 0 && simdata->params.outputs) {
      int outsampling =
          (int)(1.0 / (simdata->params.outrate * simdata->params.dt));

      printf("     Output rate: every %d step(s)\n", simdata->params.outrate);
      printf(" Output sampling: %d Hz\n\n", outsampling);
      printf(" Output files:\n\n");

      for (int i = 0; i < simdata->params.numoutputs; i++) {
        print_output(&simdata->params.outputs[i]);
      }

      printf("\n");

    } else if (simdata->params.outrate < 0) {
      printf("  Output is disabled (output rate set to 0)\n\n");

    } else {
      printf("  Output is disabled (not output specified)\n\n");
    }

    print_source(&simdata->params.source);

    fflush(stdout);
  }

  free(rho_map->vals);
  free(rho_map);
  free(c_map->vals);
  free(c_map);
}

void finalize_simulation(simulation_data_t *simdata) {
  if (simdata->params.outputs != NULL) {
    for (int i = 0; i < simdata->params.numoutputs; i++) {
      free(simdata->params.outputs[i].filename);

      if (cart_rank == root && simdata->params.outrate > 0) {
        fclose(simdata->params.outputs[i].fp);
      }
    }

    free(simdata->params.outputs);
  }

  free(simdata->params.source.data);
  free(simdata->params.cin_filename);
  free(simdata->params.rhoin_filename);

  free(simdata->rho->vals);
  free(simdata->rho);
  free(simdata->rhohalf->vals);
  free(simdata->rhohalf);
  free(simdata->c->vals);
  free(simdata->c);

  free(simdata->pold->vals);
  free(simdata->pold);
  free(simdata->pnew->vals);
  free(simdata->pnew);

  free(simdata->vxold->vals);
  free(simdata->vxold);
  free(simdata->vxnew->vals);
  free(simdata->vxnew);
  free(simdata->vyold->vals);
  free(simdata->vyold);
  free(simdata->vynew->vals);
  free(simdata->vynew);
  free(simdata->vzold->vals);
  free(simdata->vzold);
  free(simdata->vznew->vals);
  free(simdata->vznew);
  free_comm_buffers();
}

void swap_timesteps(simulation_data_t *simdata) {
  local_data_t *tmpp = simdata->pold;
  local_data_t *tmpvx = simdata->vxold;
  local_data_t *tmpvy = simdata->vyold;
  local_data_t *tmpvz = simdata->vzold;

  simdata->pold = simdata->pnew;
  simdata->pnew = tmpp;
  simdata->vxold = simdata->vxnew;
  simdata->vxnew = tmpvx;
  simdata->vyold = simdata->vynew;
  simdata->vynew = tmpvy;
  simdata->vzold = simdata->vznew;
  simdata->vznew = tmpvz;
}

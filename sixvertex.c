#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "mt19937ar.c"
#include "sixvertex.h"


#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

// From combinatorics, the maximum number of combinations of the numbers 1-6 where order 
// and repeats are not relevant is 126
#define MAX_CONFIGURATIONS 10000

// TODO: I hate global variables
int state_table[6][4][2][2] = {{{{ 3 , 1 },{ 5 , 2 },},{{ 3 , 0 },{ 4 , 3 },},{{ 5 , 0 },{ 2 , 3 },},{{ 4 , 1 },{ 2 , 2 },},},
                               {{{ 2 , 1 },{ 4 , 2 },},{{ 2 , 0 },{ 5 , 3 },},{{ 4 , 0 },{ 3 , 3 },},{{ 5 , 1 },{ 3 , 2 },},},
                               {{{ 1 , 1 },{ 5 , 3 },},{{ 1 , 0 },{ 4 , 2 },},{{ 4 , 1 },{ 0 , 3 },},{{ 5 , 0 },{ 0 , 2 },},},
                               {{{ 0 , 1 },{ 4 , 3 },},{{ 0 , 0 },{ 5 , 2 },},{{ 5 , 1 },{ 1 , 3 },},{{ 4 , 0 },{ 1 , 2 },},},
                               {{{ 1 , 2 },{ 3 , 3 },},{{ 2 , 2 },{ 0 , 3 },},{{ 1 , 0 },{ 2 , 1 },},{{ 3 , 0 },{ 0 , 1 },},},
                               {{{ 0 , 2 },{ 2 , 3 },},{{ 3 , 2 },{ 1 , 3 },},{{ 0 , 0 },{ 3 , 1 },},{{ 2 , 0 },{ 1 , 1 },},},
};

//int seed = 122;

int state_arrows[6][4] = {
    {RIGHT, RIGHT, UP, UP},
    {LEFT, LEFT, DOWN, DOWN},
    {RIGHT, RIGHT, DOWN, DOWN},
    {LEFT, LEFT, UP, UP},
    {RIGHT, LEFT, UP, DOWN},
    {LEFT, RIGHT, DOWN, UP}
};

const char direction_strs[4][10] = {"Left", "Right", "Up", "Down"};

int visited[MAX_CONFIGURATIONS] = {0};

void mark_visited(int number) {
  if (number >= 0 && number < MAX_CONFIGURATIONS) {
    visited[number] = 1;
  }
}

int is_visited(int number) {
  if (number >= 0 && number < MAX_CONFIGURATIONS) {
    return visited[number];
  }

  return 0;
}

// This function returns the offset in a 1D array
// of the coordinates of the vertex in the corresponding 2D array.
int offset2d(int x, int y, int Ny) {
  return x * Ny + y;
}

int return_vertex_state(struct vertex **vertices_ptr, int Nx, int Ny) {
  struct vertex *vertices = !vertices_ptr ? NULL : *vertices_ptr;
  if (!vertices) return -1;

  int *states = malloc(Nx*Ny*sizeof(int));
  if (!states) return -1;
  memset(states, 0, Nx*Ny*sizeof(int));

  int i = 0;
  // Go row-by-row, starting from the top
  for (int y = Ny - 1; y >= 0; y--) {
    for (int x = 0; x < Nx; x++) {
      int offset = offset2d(x, y, Ny);
      if (offset >= Nx*Ny) {
        printf("Out of bounds access into `vertices` array\n");
      }
      struct vertex *vert = vertices + offset;
      //printf("%zu\n", sizeof(struct vertex));
      int state = vert->state;
      states[i] = state;
      i++;
    }
  }
  
  int multiplier = 1;
  for (int i = 0; i < (Nx*Ny -1); i++) {
    multiplier *= 10;
  }

  int num = 0;
  int max = Nx*Ny - 1;
  //printf("%d\n", max);
  for (int i = max; i >= 0; i--) {
    //printf("index i = %d\n", i);
    //printf("%d\n", num);
    int state = states[i];
    num += (int) (multiplier * state);
    multiplier = (int)(multiplier / 10);
  }
  
  free(states);
  return num;
}

// Generate an initial lattice of dimension Nx*Ny with all vertices initially in state 2.
void init_lattice(struct vertex **vertices_ptr, int Nx, int Ny) {
  int N = Nx*Ny;

  struct vertex *vertices = !vertices_ptr ? NULL : *vertices_ptr; 
  if (!vertices) return;
  
  int xmax = Nx - 1;
  int ymax = Ny - 1;

  int c = 0;
  for (int x = 0; x < Nx; x++) {
    for (int y = 0; y < Ny; y++) {
      struct vertex *v = vertices + c;
      // Determine if this is a special vertex (on left or right edge) or not.
      // 1 -> left edge
      // 2 -> right edge
      // 0 -> neither
      if (!x) v->pos_x = 1;
      else if (x == xmax) v->pos_x = 2;
      else v->pos_x = 0;
      // Determine if this is a special vertex (on top or bottom edge).
      // Notation is analogous to the x case
      if (!y) v->pos_y = 1;
      else if (y == ymax) v->pos_y = 2;
      else v->pos_y = 0;
      v->row = y;
      v->col = x;

      v->state = 1;     
      if (x > 0) v->neighbours[LEFT] = vertices + offset2d(x - 1, y, Ny);
      else v->neighbours[LEFT] = NULL;
      
      if (x < xmax) v->neighbours[RIGHT] = vertices + offset2d(x + 1, y, Ny);
      else v->neighbours[RIGHT] = NULL;

      if (y > 0) v->neighbours[DOWN] = vertices + offset2d(x, y - 1, Ny);
      else v->neighbours[DOWN] = NULL;

      if (y < ymax) v->neighbours[UP] = vertices + offset2d(x, y + 1, Ny);
      else v->neighbours[UP] = NULL;
      
      c++;
    }
  }

  return;
}

// Read a particular state of the lattice from file.
void read_lattice_from_file(struct vertex **vertices_ptr, char *filename, int Nx, int Ny) {
  int N = Nx * Ny;
  struct vertex *vertices = !vertices_ptr ? NULL : *vertices_ptr; 
  
  printf("%s\n", filename);
  FILE *fptr = fopen(filename, "r");
  
  char *buf = NULL;
  size_t bufsize = 0;
   
  int res = getline(&buf, &bufsize, fptr);
  if (res == -1) {
    printf("Error: Malformed data file for initial lattice state. Exiting..\n");
    exit(1);
  }

  int xmax = Nx - 1;
  int ymax = Ny - 1;

  while (getline(&buf, &bufsize, fptr) != -1) {
    int x, y, state, idx;
    int tot = sscanf(buf, "%d %d %d %d", &x, &y, &state, &idx);
    if (tot < 4) {
      printf("Error reading lattice positions from file %s", filename);
    }
    
    struct vertex *v = vertices + idx;
    //printf("The vertex at position %d has coordinates (%d, %d) and is in state %d\n", idx, x, y, state);

    // Determine if this is a special vertex (on left or right edge) or not.
    // 1 -> left edge
    // 2 -> right edge
    // 0 -> neither
    if (!x) v->pos_x = 1;
    else if (x == xmax) v->pos_x = 2;
    else v->pos_x = 0;
    // Determine if this is a special vertex (on top or bottom edge).
    // Notation is analogous to the x case
    if (!y) v->pos_y = 1;
    else if (y == ymax) v->pos_y = 2;
    else v->pos_y = 0;
    v->row = y;
    v->col = x;

    v->state = state;

    if (x > 0) v->neighbours[LEFT] = vertices + offset2d(x - 1, y, Ny);
    else v->neighbours[LEFT] = NULL;
      
    if (x < xmax) v->neighbours[RIGHT] = vertices + offset2d(x + 1, y, Ny);
    else v->neighbours[RIGHT] = NULL;

    if (y > 0) v->neighbours[DOWN] = vertices + offset2d(x, y - 1, Ny);
    else v->neighbours[DOWN] = NULL;

    if (y < ymax) v->neighbours[UP] = vertices + offset2d(x, y + 1, Ny);
    else v->neighbours[UP] = NULL;
  }

  if (feof(fptr)) {
    printf("Finished reading vertex state.\n");
  }
   
  //check_system(vertices, N);
  fclose(fptr);
  free(buf);

  return;
}

// Calculate the total energy of the system according to whether the spins
// at the rightmost edge are aligned parallel or anti-parallel with the 
// direction of the rightmost magnet, whose alignment is determined by 
// the sign of J. This function calculates the ratio of the total interface 
// energy to the thermal energy kBT by multiplying by the dimensionless ratio 
// `betaJ`, where `beta` is thermodynamic beta 1/kBT
double calculate_total_energy(struct vertex *vertices, double betaJ, int N) {
  double E = 0.0;
  if (!vertices) {
    printf("Error: Got NULL pointer in calculate_total_energy(). Exiting program...\n");
    exit(0);
  }

  for (int i = 0; i < N; i++) {
    struct vertex *vert = vertices + i;
    
    if (vert->pos_x == 2) {
      int align = state_arrows[vert->state][RIGHT] == RIGHT;
      // By convention, treat RIGHT as +1 and LEFT as -1.
      E += betaJ * (2*align - 1); 
    }
  }

  return E / N;
}

// This is the quantity that we are interested in;
// integrated over values of the coupling it gives the free energy 
// of the system.
double calculate_interface_magnetisation(struct vertex *vertices, int N) {
  double m = 0.0;
  
  for (int i = 0; i < N; i++) {
    struct vertex *vert = vertices + i;
    if (vert->pos_x == 2) {
      int align = state_arrows[vert->state][RIGHT] == RIGHT;
      int sign = 2*align - 1; // = +/- 1
      m += sign;
    }
  }

  return m / N;
}

// The vertex density is a good metric for the equilibration of the 
// six vertex model - see Newman & Barkema: Monte Carlo Methods in Statistical Physics
double calculate_vertex_density(struct vertex *vertices, int N) {
  int rho = 0;

  for (int i = 0; i < N; i++) {
    struct vertex *vert = vertices + i;
    if (vert->state == 4 || vert->state == 5) rho += 1;
  }

  return rho / (double)N;
}



struct vertex *move(struct vertex *vertices, double betaJ, int N, int Nx, int Ny, int *state_count) {
  struct vertex *backup_vertices = calloc(N, sizeof(struct vertex));
  memcpy(backup_vertices, vertices, N * sizeof(struct vertex));

  // the opposites of the four arrow directions, stored in the right order based
  // on the definitions of UP, DOWN, LEFT and RIGHT such that e.g. opposite[RIGHT] = LEFT.
  int opposite[4] = {RIGHT, LEFT, DOWN, UP};
  
  // Randomly pick a vertex to begin this move from.
  int index = genrand_real2() * N;
  
  // Preserve the initial state of the initially chosen vertex.
  struct vertex *initial_vertex = vertices + index;
  int initial_state = initial_vertex->state;
  
  // Randomly generate a direction to move in 
  // and an arrow to flip.
  int dir = genrand_real2() * 4;
  int choice = genrand_real2() * 2;
  
  // The new state of the chosen vertex after this choice of spin flip.
  int new_state = state_table[initial_state][dir][choice][0];
  // The other direction that should be flipped to preserve the ice rules.
  int opp_dir = state_table[initial_state][dir][choice][1];
  // The alternative state of the chosen vertex had we made the opposite choice.
  int alt_state = state_table[initial_state][dir][1 - choice][0];

  initial_vertex->state = new_state;
  
  // Store the current vertex, previous vertex, and current direction the move is taking.
  struct vertex *current = initial_vertex;
  int curr_dir = dir;
  int done = 0;
  struct vertex *last = initial_vertex;

  while (!done) {
    last = current;
    // move in the specified direction.
    current = current->neighbours[curr_dir];
    if (!current) {
      // the vertex does not have a neighbour in the direction we are moving in, corresponding to hitting an edge.
      done = 2;
      break;
    }
    if (current == initial_vertex) {
      done = 1;
      break;
    }
    
    // Randomly choose where to move next.
    choice = genrand_real2() * 2;
    curr_dir = opposite[curr_dir];
    int *new_state_data = state_table[current->state][curr_dir][choice];
    current->state = new_state_data[0];
    curr_dir = new_state_data[1];
  } 
  
  if (done == 1) {
    int in_dir = opposite[curr_dir];
    if (in_dir == dir) initial_vertex->state = initial_state;
    else if (in_dir != opp_dir) initial_vertex->state = alt_state;
    free(backup_vertices);
    return vertices;
  }

  if (curr_dir == LEFT) {
    // We hit an edge, and we are moving to the left. We must be at the left of the lattice.
    memcpy(vertices, backup_vertices, N*sizeof(struct vertex));
    free(backup_vertices);
    return vertices;
  }
  
  double delta_E = 0;

  if (curr_dir == RIGHT) {
    if (state_arrows[last->state][RIGHT] == RIGHT) delta_E += -2*betaJ;
    else delta_E += 2*betaJ;
  } 
  
  // We have not managed to close the loop yet. Return to the initial vertex, and move in the opposite direction instead.
  // Do this until we hit another edge, which is then a valid move.
  current = initial_vertex;
  curr_dir = opp_dir;
  done = 0;
  
  // Go back along the opposite direction until we hit another edge.
  while (!done) {
    last = current;
    current = current->neighbours[curr_dir];
    if (!current) {
      done = 2;
      break;
    }
    // Again, randomly decide which direction to move in.
    choice = genrand_real2() * 2;
    curr_dir = opposite[curr_dir];
    int *new_state_data = state_table[current->state][curr_dir][choice];
    current->state = new_state_data[0];
    curr_dir = new_state_data[1];
  } 

  if (curr_dir == LEFT) {
    // We hit an edge, and we are moving to the left. We must be at the left of the lattice.
    //printf("I hit the left edge. %d %d", last->pos_x, last->pos_y);
    memcpy(vertices, backup_vertices, N*sizeof(struct vertex));
    free(backup_vertices);
    return vertices;
  }

  if (curr_dir == RIGHT) {
    if (state_arrows[last->state][RIGHT] == RIGHT) delta_E += -2*betaJ;
    else delta_E += 2*betaJ;
  }
   
  if (genrand_real2() > exp(-delta_E)) {
    memcpy(vertices, backup_vertices, N * sizeof(struct vertex));
    free(backup_vertices);
    // reject the move 
    return vertices; 
  }
  
  int state = return_vertex_state(&vertices, Nx, Ny); 
  if (!is_visited(state)) {
    *state_count += 1;
    mark_visited(state);
  }
  
  free(backup_vertices);
  // Accept the move.
  return vertices;
}

void write_lattice_state(struct vertex *vertices, char *filename, int Nx, int Ny) {
  int N = Nx * Ny;

  FILE *fptr = fopen(filename, "w");
  if (!fptr) return;
  
  printf("Writing lattice state to file %s\n", filename);
  fprintf(fptr, "x y state idx\n");
  
  for (int i = 0; i < N; i++) {
    struct vertex *vert = vertices + i;
    int x = vert->col;
    int y = vert->row;
    int state = vert->state;

    fprintf(fptr, "%d %d %d %d\n", x, y, state, i);
  } 

  fclose(fptr);
}

void check_vertex_consistency(struct vertex *v) {
  int opposite[4] = {RIGHT, LEFT, DOWN, UP};

  for (int i = 0; i < 4; i++) {
    if (!v->neighbours[i]) continue;
    if (state_arrows[v->neighbours[i]->state][opposite[i]] != state_arrows[v->state][i]) {
      printf("Error! State %d to the %s of %d at position (%d, %d)\n", v->neighbours[i]->state, direction_strs[i], v->state, v->col, v->row);
      exit(3);
    }
  }
}

void check_system(struct vertex *vertices, int N) {
  for (int i = 0; i < N; i++) {
    check_vertex_consistency(vertices + i);
  }
}

int main(int argc, char *argv[]) {
  //int r  = rand();
  //init_genrand(r);
  
  if (argc < 7) {
    printf("six_vertex: Requires 6 arguments, but %d were given\n", (argc - 1));
    printf("Argument format: ./sixvertex betaJ N_x N_y outfile equ_file \n");
    printf("betaJ: Ratio of interfacial exchange coupling to thermal energy.\n");
    printf("N_x: Number of sites in x direction.\n");
    printf("N_y: Number of sites in y direction.\n");
    printf("outfile: File to write output to.\n");
    printf("equ_file: File that the final (equilibrated) state of the system will be written to.\n");
    printf("seed_file: File that the RNG seed used by this simulation is \n");
    return -1;
  }
  
  // Read in parameters as command line arguments.
  double betaJ = strtod(argv[1], NULL);
  int Nx = strtol(argv[2], NULL, 10);
  int Ny = strtol(argv[3], NULL, 10);
  int N = Nx * Ny;
  
  // Initialise random number generator.
  srand(time(0));
  int r = rand();
  init_genrand(r); 

  char *outfile = argv[4];
  char *equ_file = argv[5]; 
  char *seed_file = argv[6];
 
  FILE *seed_fptr = fopen(seed_file, "w");
  fprintf(seed_fptr, "SEED USED:\n");
  fprintf(seed_fptr, "%d", r);
  fclose(seed_fptr);
  
  printf("six_vertex: Using below parameters\n");
  printf("betaJ: %f\n", betaJ);
  if (betaJ >= 0) {
    printf("J > 0\n");
  } else {
    printf("J < 0\n");
  }
  printf("N_x: %d\n", Nx);
  printf("N_y: %d\n", Ny);
  printf("outfile: %s\n", outfile);
  printf("equ_file: %s\n", equ_file);
  
  struct vertex *vertices = calloc(N,sizeof(struct vertex));
  printf("size of vertex data in KB %f\n", (N*sizeof(struct vertex)) / ((double)1024));
  
  init_lattice(&vertices, Nx, Ny);

  //printf("This output occurs after reading the data from file\n");
  //write_lattice_state(vertices, "test_data.vertex", Nx, Ny);
  check_system(vertices, N);
  
  int calc_steps = 5e5;
  int equ_steps = 1e7;
  int output_steps = 20;
  int equ_output_steps = 1;
  int flush_steps = 1e3;
   
  FILE *fptr = fopen(outfile, "w");
  fprintf(fptr, "iter  energy  m_int  rho\n");
  
  clock_t start, end;
  double cpu_time;
  
  // Equilibrate the system.
  printf("Equilibrating for %e steps...\n", (double)equ_steps);
  int state_count = 0;
  for (int i = 0; i < equ_steps; i++) {
    vertices = move(vertices, betaJ, N, Nx, Ny, &state_count);
   // if (!(i % equ_output_steps)) {
   //   double E = calculate_total_energy(vertices, betaJ, N);
   //   double m = calculate_interface_magnetisation(vertices, N);
   //   double rho = calculate_vertex_density(vertices, N);
   //   fprintf(fptr, "%d  %f  %f  %f\n", i, E, m, rho);
   //   fflush(fptr);
   // }
  }

  //printf("This output occurs after equilibrating the system\n");
  write_lattice_state(vertices, equ_file, Nx, Ny);
  check_system(vertices, N);

  printf("Calculating thermodynamic properties for %e steps...\n", (double)calc_steps);
  for (int i = 0; i < calc_steps; i++) {
    vertices = move(vertices, betaJ, N, Nx, Ny, &state_count);
    if (!(i % output_steps)) {
      double E = calculate_total_energy(vertices, betaJ, N);
      double m = calculate_interface_magnetisation(vertices, N);
      double rho = calculate_vertex_density(vertices, N);
      fprintf(fptr, "%d  %f  %f  %f\n", i, E, m, rho);
    }
  }
  fclose(fptr);
  free(vertices);

  return 0;
}

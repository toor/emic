#ifndef __SIXVERTEX_H__
#define __SIXVERTEX_H__

struct vertex {
  int state;
  struct vertex *neighbours[4];
  int pos_x;
  int pos_y;
  int row;
  int col;
};

struct vertex *move(struct vertex *, double, int, int, int, int *);

double calculate_total_energy(struct vertex *, double, int);
double calculate_total_magnetisation(struct vertex *, int, int);

void write_lattice_state(struct vertex *, char *, int, int);

int offset2d(int, int, int);

void init_lattice(struct vertex **, int, int);
void read_lattice_from_file(struct vertex **, char *, int, int);

void check_vertex_consistency(struct vertex *);
void check_system(struct vertex *, int);

#endif

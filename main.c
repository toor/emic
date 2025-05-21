#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

#define L 3

#define N (L*L)
#define UP 2
#define DOWN (-UP)
#define RIGHT 3
#define LEFT (-RIGHT)

struct vertex {
  int arrows[4];
  int type;
};


double rand01(void) {
  return (double)rand() / (double)((unsigned)RAND_MAX + 1);
} 

bool is_empty_vertex(struct vertex *vertex) {
  if (vertex->type == -1) return true;
  else return false;
}

int rand_range(int limit) {
  int divisor = RAND_MAX / (limit + 1);
  int retval;

  do {
    retval = rand() / divisor;
  } while (retval > limit);

  return retval;
}

int offset3d(int i, int j, int k) {
  return (k * L * L) + (j * L) + i;
}


/* Determines whether we have visited a particular site before in this iteration 
 * Params:
 * int *steps: A list of steps that have been taken in the current iteration.
 * int *lasttime: Stores the  
 *
 * */
int last_visited(int **steps, int i, int j, int k, int len) {
  int m = len - 1;
  /* Count down from the current iteration */
  for (; m >= 0; m--) {
    if (steps[m][0] == i && steps[m][1] == j && steps[m][2] == k) {
      return m;
    }
  }
}

void reverse_arrows(struct vertex **vertices, int **steps, int start, int len) {
  int m = len -1;
  for (; m >= 0; m--) {
    int arrow = steps[m][2];
    int i = steps[m][0];
    int j = steps[m][1];

    switch (arrow) {
      case 0:
      // Arrow pointing up; reverse the down arrow on the vertex above as well as the up arrow on this vertex.
      vertices[i][j].arrows[arrow] *= -1;
      vertices[i - 1][j].arrows[arrow] *= -1;
      break;
      case 1:
      // Arrow pointing right; reverse the left arrow on the adjacent arrow.
      vertices[i][j].arrows[arrow] *= -1;
      vertices[i][j+1].arrows[arrow] *= -1;
      break;
      case 2:
      // Arrow pointing down; reverse the up arrow on the adjacent vertex.
      vertices[i][j].arrows[arrow] *= -1;
      vertices[i+1][j].arrows[arrow] *= -1;
      break;
      case 3:
      // Arrow pointing left; reverse the left arrow on the adjacent vertex.
      vertices[i][j].arrows[arrow] *= -1;
      vertices[i][j-1].arrows[arrow] *= -1;
    }
  }
}

/* This function actually performs the Monte Carlo move each step. 
* Each step, we need to determine which type of vertex we have; This
* also determines which sides of the vertex the outgoing arrows come from.
* We then store a list of the coordinates i,j of the chosen vertex, which sides
* the reversal takes place on, and which side the "fixing" of this defect takes place on.
* EXAMPLE: A vertex of type 1 has outgoing arrows up and to the right. If we randomly choose
* to reverse the arrow pointing up, we can either reverse the arrow at the bottom
* or the arrow at the left. If we choose to reverse the arrow at the left,
* we store the array [i, j, 1, 4] as the next element in the list of steps.
* The different arrows are labelled as follows:
* Above => State 1
* Right => State 2
* Below => State 3
* Left => State 4
*/
void monte_carlo_move(struct vertex **vertices, int **steps, int **lasttime, int iter, double beta, double J) {
  int i, i0 = rand_range(L - 1);
  int j, j0 = rand_range(L - 1);
  int start_arrow;

  int len = 0, noloop = 1;
  // This variable should be set if we hit the rightmost edge of the lattice, and therefore modify the energy
  // of interaction with the right magnet.
  int is_edge = 0;

  do {
    // j is the column index; if we hit the leftmost column we should reject the move.
    if (!j) {
      return;
    }

    if (i == 0 || i == (L - 1) || j == (L - 1)) {
      if (!len) continue; // skip the trivial case
      if (i0 == 0 || i0 == L - 1 || j0 == (L - 1)) {
        noloop = 0; // end the loop; we have found a valid configuration
        is_edge = 1;
        continue;
      } else {
        return; // we hit the edge, but didn't start at another edge. This is not a valid move.
      }
    } else if (lasttime[i][j] == iter) {
      noloop = 0;
      continue;
    }
    
    switch (vertices[i][j].type) {
      case 1:
        // Outgoing arrows are to the right and up
        if (rand01() < 0.5) {
          start_arrow = 1;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j += 1; // move to the right.
          continue;
        } else {
          start_arrow = 0;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i -= 1; // move up
          continue;
        }
        break;
      case 2:
        // Outgoing arrows are to the left and down.
        if (rand01() < 0.5) {
          start_arrow = 2;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i += 1; // move down.
          continue;
        } else {
          start_arrow = 3;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j -= 1; // move left.
          continue;
        }
        break;
      case 3:
        // Outgoing arrows are to the right and down.
        if (rand01() < 0.5) {
          start_arrow = 1;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j += 1; // move right.
          continue;
        } else {
          start_arrow = 2;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i += 1;
          continue;
        }
        break;
      case 4:
        // Outgoing arrows are to the left and up.
        if (rand01() < 0.5) {
          start_arrow = 0;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i -= 1; // move up.
          continue;
        } else {
          start_arrow = 3;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j -= 1; // move left.
          continue;
        }
        break;
      case 5:
        // Outgoing arrows are up and down
        if (rand01() < 0.5) {
          start_arrow = 0;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i -= 1;
          continue;
        } else {
          start_arrow = 2;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          i += 1;
          continue;
        }
        break;
      case 6:
        // Outgoing arrows are left and right.
        if (rand01() < 0.5) {
          start_arrow = 1;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j += 1; // move to the right.
          continue;
        } else {
          start_arrow = 3;
          steps[len++][0] = i;
          steps[len++][1] = j;
          steps[len++][2] = start_arrow;
          lasttime[i][j] = iter;
          j -= 1; // move to left.
          continue;
        }
        break;
    }
  } while (noloop);

  if (is_edge) {
    if (j == L - 1) {
      double E_0, E_1 = 0.0;

      for (int m = 0; m < L; m++) {
        E_0 += vertices[L-1][m].arrows[1];
      }
      
      // reverse the arrow to the right as well.
      vertices[i][j].arrows[1] *= -1;
      reverse_arrows(vertices, steps, 0, len);


      for (int m = 0; m < L; m++) {
        E_1 += vertices[L-1][m].arrows[1];
      }
      // TODO accept on Boltzmann distribution.
      double delta_E = J * (E_1 - E_0);
      if (delta_E < 0.0) {
        // accept move
        goto write_output;
      } else if (rand01() < exp(-beta * delta_E)) {
        // accept move
        goto write_output;
      } else {
        vertices[i][j].arrows[1] *= -1;
        reverse_arrows(vertices, steps, 0, len);
      }
    write_output:
      if (iter % output_steps) {
        FILE *fptr = fopen(output_file, "w");
        fprintf(fptr, "%d", iter);
        fprintf(fptr, "")
        fprintf(fptr, "%f", E_1);
        fprintf(fptr, "\n");
        fclose(fptr);
      }
      return;
    } else if (i == 0) {
      // reverse the top arrow as well.
      vertices[i][j].arrows[0] *= -1;
      reverse_arrows(vertices, steps, 0, len); 
    } else if (i == L - 1) {
      // reverse the bottom arrow as well.
      vertices[i][j].arrows[2] *= -1;
      reverse_arrows(vertices, steps, 0, len); 
    }
  } else {
    // lasttime[i][j] has been set during the current Monte Carlo step. Work backwards through the list of arrow moves
    // to determine the step where it was last visited (i.e. check whether a short or long loop has been formed).
    int m = last_visited(steps, i, j, start_arrow, len);
    
    reverse_arrows(vertices, steps, m, len);
  }
}


void fill_arrows(struct vertex *vert, int arr[4]) {
  for (int i = 0; i < 4; i++) {
    vert->arrows[i] = arr[i];
  }
}

void generate_initial_lattice_state(struct vertex **vertices, char *outfile) { 
  int TYPE1[4] = {UP, RIGHT, UP, RIGHT};
  int TYPE2[4] = {DOWN, LEFT, DOWN, LEFT};
  int TYPE3[4] = {DOWN, RIGHT, DOWN, RIGHT};
  int TYPE4[4] = {RIGHT, UP, RIGHT, UP};
  int TYPE5[4] = {UP, LEFT, DOWN, RIGHT};
  int TYPE6[4] = {DOWN, RIGHT, UP, LEFT};


  for (int i = 0; i < L; i++) {
    printf("Filling in first column.");
    vertices[i][0].type = 1;
    fill_arrows(&vertices[i][0], TYPE1);
  }
  
  // fill in first row randomly
  for (int i = 1; i < L; i++) {
    int prev = vertices[0][i - 1].type;
    
    // arrow pointing to the right
    if (prev == 1 || prev == 3 || prev == 6) {
      // Choose between 1, 2 or 3 randomly.
      // 1 -> 1, 2 -> 3, 3 -> 5.
      int choice = rand() % (3 - 1 + 1) + 1;
      switch (choice) {
        case 1:
          vertices[0][i].type = 1;
          fill_arrows(&vertices[0][i], TYPE1);
          break;
        case 2:
          vertices[0][i].type = 3;
          fill_arrows(&vertices[0][i], TYPE3);
          break;
        case 3:
          vertices[0][i].type = 5;
          fill_arrows(&vertices[0][i], TYPE5);
          break;
      }
    } else if (prev == 2 || prev == 4 || prev == 5) {
      // 1 -> 2, 2 -> 4, 3 -> 6
      int choice = rand() % (3 - 1 + 1) + 1;
      switch (choice) {
        case 1:
          vertices[0][i].type = 2;
          fill_arrows(&vertices[0][i], TYPE2);
          break;
        case 2:
          vertices[0][i].type = 4;
          fill_arrows(&vertices[0][i], TYPE4);
          break;
        case 3:
          vertices[0][i].type = 6;
          fill_arrows(&vertices[0][i], TYPE6);
          break;
      }
    }
  }

  for (int i = 1; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int vert_left = vertices[i][j-1].type;
      int vert_above = vertices[i-1][j].type;

      if (vert_left == 1 || vert_left == 3 || vert_left == 6) {
        if (vert_above == 1 || vert_above == 4 || vert_above == 6) {
          if (rand01() < 0.5) {
            vertices[i][j].type = 1;
            fill_arrows(&vertices[i][j], TYPE1);
          } else {
            vertices[i][j].type = 5;
            fill_arrows(&vertices[i][j], TYPE5);
          } 
        } else {
          vertices[i][j].type = 3;
          fill_arrows(&vertices[i][j], TYPE3);
        }
      } else if (vert_left == 2 || vert_left == 4 || vert_left == 5) {
        if (vert_above == 2 || vert_above == 3 || vert_above == 5) {
          if (rand01() < 0.5) {
            vertices[i][j].type = 2;
            fill_arrows(&vertices[i][j], TYPE2);
          } else {
            vertices[i][j].type = 6;
            fill_arrows(&vertices[i][j], TYPE6);
          }
        } else {
          vertices[i][j].type = 4;
          fill_arrows(&vertices[i][j], TYPE4);
        }   
      }
    }
  }

  FILE *fptr = fopen(outfile, "w");

  char *header = "row,col,up,right,down,left\n";
  fprintf(fptr, "%s", header);

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      struct vertex vert = vertices[i][j];

      fprintf(fptr, "%d", i);
      fprintf(fptr, ",");
      fprintf(fptr, "%d", j);
      fprintf(fptr, ",");
      fprintf(fptr, "%d", vert.arrows[0]);
      fprintf(fptr, ",");
      fprintf(fptr, "%d", vert.arrows[1]);
      fprintf(fptr, ",");
      fprintf(fptr, "%d", vert.arrows[2]);
      fprintf(fptr, ",");
      fprintf(fptr, "%d", vert.arrows[3]);
      fprintf(fptr, "\n");
    }
  }
}


int main(int argc, char *argv[]) {
  if (!argc) {
    printf("Error: This program takes two arguments: energetic multiplier and temperature. However, no arguments were provided. Exiting...");
    return -1;
  } 

  double J = strtod(argv[0], NULL);
  double T = strtod(argv[1], NULL);
  double beta = 1/T;
  int MAX_STEPS = 1000000;
  int output_steps = 1000;

  struct vertex **vertices = calloc(L, sizeof(struct vertex *));
  int **lasttime = calloc(L, sizeof(int *));

  if (vertices == NULL) {
    printf("Error allocating requisite memory for this program. Exiting...");
    return -1;
  }

  for (int i = 0; i < L; i++) {
    vertices[i] = calloc(L, sizeof(struct vertex));
    if (vertices[i] == NULL) {
      printf("Error allocating requisite memory for this program. Exiting...");
      return -1;
    }
    lasttime[i] = calloc(L, sizeof(int));
    if (lasttime[i] == NULL) {
      printf("Error allocating requisite memory for this program. Exiting...");
    }
  }
   
  char *outfile = "initial_lattice.csv";
  generate_initial_lattice_state(vertices, outfile);
  

  int **steps = calloc(N, sizeof(int *));
  if (steps == NULL) {
    printf("Error allocating requisite memory for this program. Exiting...");
    return -1;
  }

  for (int i = 0; i < N; i++) {
    steps[i] = calloc(3, sizeof(int));
    if (steps[i] == NULL) {
      printf("Error allocating requisite memory for this program. Exiting...");
      return -1;
    }
  }
  
  printf("Begin equilibration period at a temperature of %f", T);
  for (int i = 0; i < MAX_STEPS; i++) {
    monte_carlo_move(vertices, steps, lasttime, i, J, beta);
  }
  
  return 0;
}
//
//void move(int *arrow, int *step, int *lasttime, int iter) {
//  int i;
//  int first, curr;
//  int currarrow;
//  int nc,cand;
//  int *clist = malloc(2*sizeof(int));
//
//
//  int len=0,noloop=1;
//
//  first = N*rand();
//  lasttime[first] = iter;
//  step[len++] = curr = first;
//  currarrow = arrow[curr];
//
//  // Generate loop
//  do {
//    nc = 0;
//
//    switch (currarrow) {
//      case LU:
//        if ((cand=curr-XNN)<0) cand += N;
//        if (arrow[cand]==LD) clist[nc++] = cand;
//        if ((cand=curr-(XNN+YNN)) < 0) cand+= N;
//        if (arrow[cand]==LU) clist[nc++] = cand;
//        if ((cand=curr-YNN) < 0) cand += N;
//        if (arrow[cand]==RU) clist[nc++] = cand;
//        break;
//      case LD:
//        if ((cand=curr-XNN) < 0) cand += N;
//        if (arrow[cand]==LU) clist[nc++] = cand;
//        if ((cand=curr-XNN+YNN)>=N) cand -= N;
//        if (arrow[cand]==LD) clist[nc++] = cand;
//        if ((cand=curr+YNN)>=N) cand -= N;
//        if (arrow[cand]==RD) clist[nc++] = cand;
//        break;
//      case RU:
//        if ((cand=curr+XNN)>=N) cand -= N;
//        if (arrow[cand]==RD) clist[nc++] = cand;
//        if ((cand=curr+XNN-YNN)<0) cand += N;
//        if (arrow[cand]==RU) clist[nc++] = cand;
//        if ((cand=curr-YNN)<0) cand += N;
//        if (arrow[cand]==LU) clist[nc++] = cand;
//        break;
//      case RD:
//        if ((cand=curr+XNN)>=N) cand -= N;
//        if (arrow[cand]==RU) clist[nc++] = cand;
//        if ((cand=curr+XNN+YNN)>=N) cand -= N;
//        if (arrow[cand]==RD) clist[nc++] = cand;
//        if ((cand=curr+YNN)>=N) cand -= N;
//        if (arrow[cand]==LD) clist[nc++] = cand;
//        break;
//    }
//
//    if (lasttime[clist[0]] == iter) {
//      step[len++] = clist[0];
//      noloop = 0;
//    } else if (lasttime[clist[1]] == iter) {
//      step[len++] = clist[1];
//      noloop = 0;
//    } else {
//      if ((float)rand() < 0.5) curr = clist[0];
//      else curr = clist[1];
//      lasttime[curr] = iter;
//      step[len++] = curr;
//      currarrow = arrow[curr];
//    }
//  } while (noloop);
//
//  // Follow the loop and reverse all arrows in iter
//  for (i = 0; step[i] != step[len - 1]; i++);
//  for (; i<len-1; i++) arrow[step[i]] *= -1;
//}
//
//
//// Generate an initial lattice state.
//void generate_initial_lattice_state(int **arrows, struct vertex **vertices) {
//  char *arrows_filename = "arrows.csv";
//  char *vertices_filename = "vertices.csv";
//  
//  FILE *arrow_file;
//  FILE *vertex_file;
//
//  struct vertex empty_vert = {{-1, -1}, {-1, -1}, -1};
//
//  int initial_arrow = (rand() % (3 - 2 + 1)) + 2;
//  arrows[0][0] = initial_arrow;
//  int prev_arrow = initial_arrow;
//  
//  // Generate the first row of arrows
//  for (int i = 1; i < L; i++) {
//    if (prev_arrow == LU || prev_arrow == RD) {
//      if (rand01() < 0.5) arrows[0][i] = LD;
//      else arrows[0][i] = RU;
//    } else if (prev_arrow == RU || prev_arrow == LD) {
//      if (rand01() < 0.5) arrows[0][i] = RD;
//      else arrows[0][i] = LU;
//    }
//    prev_arrow = arrows[0][i];
//  }
//
//  // Generate the first row of vertices
//  for (int i = 0; i < L; i++) {
//    // Only fill in even-numbered vertices.
//    if (i % 2) {
//      struct vertex vert = {{arrows[0][i], arrows[0][i+1]}, {-1, -1}, -1};
//
//      vertices[0][i] = vert;
//    } else {
//      vertices[0][i] = empty_vert;
//    }
//  }
//
//  // Rows
//  for (int i = 1; i < L; i++) {
//    // On an odd (previous even) row; fill in only odd-numbered vertices
//    if ((i - 1) % 2) {
//      for (int j = 0; j < L; j += 2) {
//        struct vertex vert_above = vertices[i - 1][j];
//        struct arrows above_arrows = vert_above.above;
//        
//        struct arrows below_arrows = {-1,-1};
//        
//        int type;
//
//        // fill in the arrows below the vertex, respecting the ice rules.
//        if (above_arrows.left == RD && above_arrows.right == RU) {
//          // cases 1 and 5
//          if (rand01() < 0.5) {
//            below_arrows.left = arrows[i][j]= RU;
//            below_arrows.right = arrows[i][j+1] = RD;
//            type = 1;
//          } else {
//            below_arrows.left = arrows[i][j] = LD;  
//            below_arrows.right = arrows[i][j+1] = LU;
//            type = 5;
//          }
//        } else if (above_arrows.left == LU && above_arrows.right == LD) {
//          // cases 2 and 6
//          if (rand01() < 0.5) {
//            below_arrows.left =arrows[i][j] = LD;
//            below_arrows.right = arrows[i][j+1] = LU;
//            type = 2;
//          } else {
//            below_arrows.left = arrows[i][j] = RU;
//            below_arrows.right = arrows[i][j+1] = RD;
//            type = 6;
//          }
//        } else if (above_arrows.left==RD && above_arrows.right==LD) {
//          // case 3
//          below_arrows.left = arrows[i][j] = LD;
//          below_arrows.right = arrows[i][j+1] = RD;
//          type = 3;
//        } else if (above_arrows.left == LU && above_arrows.right == RU) {
//          // case 4
//          below_arrows.left = arrows[i][j] = RU;
//          below_arrows.right = arrows[i][j+1] = LU;
//          type = 4;
//        }
//        vert_above.below = below_arrows;
//        vert_above.type = type;
//      }
//    } else {
//      // On an odd row, only fill in odd-numbered vertices.
//      for (int j = 1; j < L; j += 2) {
//        struct vertex vert_above = vertices[i - 1][j];
//        struct arrows above_arrows = vert_above.above;
//        
//        struct arrows below_arrows = {-1,-1};
//
//        int type;
//
//        // fill in the arrows below the vertex.
//        if (above_arrows.left == RD && above_arrows.right == RU) {
//          // cases 1 and 5
//          if (rand01() < 0.5) {
//            below_arrows.left = arrows[i][j] = RU;
//            below_arrows.right = arrows[i][j+1] = RD;
//            type = 1;
//          } else {
//            below_arrows.left = arrows[i][j] = LD;
//            below_arrows.right = arrows[i][j+1] = LU;
//            type = 5;
//          }
//        } else if (above_arrows.left == LU && above_arrows.right == LD) {
//          // cases 2 and 6
//          if (rand01() < 0.5) {
//            below_arrows.left = arrows[i][j] = LD;
//            below_arrows.right = arrows[i][j+1] = LU;
//            type = 2;
//          } else {
//            below_arrows.left = arrows[i][j] = RU;
//            below_arrows.right = arrows[i][j+1] = RD;
//            type = 6;
//          }
//        } else if (above_arrows.left==RD && above_arrows.right==LD) {
//          // case 3
//          below_arrows.left = arrows[i][j] = LD;
//          below_arrows.right = arrows[i][j+1] = RD;
//          type = 3;
//        } else if (above_arrows.left == LU && above_arrows.right == RU) {
//          // case 4
//          below_arrows.left = arrows[i][j+1] = RU;
//          below_arrows.right = arrows[i][j+1] = LU;
//          type = 4;
//        }
//        vert_above.below = below_arrows;
//        vert_above.type = type;
//      }
//    }
//  }
//
//  arrow_file = fopen(arrows_filename, "w");
//
//  for (int i = 0; i < L; i++) {
//    for (int j = 0; j < L; j++) {
//      fprintf(arrow_file, "%d", arrows[i][j]);
//      fprintf(arrow_file, ",");
//    }
//    fprintf(arrow_file, "\n");
//  }
//  fclose(arrow_file);
//
//  vertex_file = fopen(vertices_filename, "w");
//  for (int i = 0; i < L; i++) {
//    for (int j = 0; j < L; j++) {
//      struct vertex vert = vertices[i][j];
//      if (is_empty_vertex(&vert)) {
//        fprintf(vertex_file, "%d", 0); // 0 encodes an empty vertex
//        fprintf(vertex_file, ",");
//      } else {
//        fprintf(vertex_file, "%d", vert.type);
//        fprintf(vertex_file, ",");
//      }
//    }
//    fprintf(vertex_file, "\n");
//  }
//
//  fclose(vertex_file);
//}
//
//
//// Write the entire state of the rotated lattice to a file.
//int store_lattice_state(int *arrows) {
//  FILE *file = fopen("output.csv", "w");
//
//  if (file == NULL) return -1;
//
//  for (int i = L; i < N; i += L) {
//    int *row = malloc(L*sizeof(int));
//    if (row == NULL) return -1;
//
//    for (int j = 0; j < i; j++) {
//      row[j] = arrows[j];
//    }
//
//    for (int k = 0; k < L; k++) {
//      fprintf(file, "%d", row[k]);
//
//      if (k < L - 1) {
//        fprintf(file, ",");
//      }
//    }
//    fprintf(file, "\n");
//    free(row);
//  }
//}
//
//int main(int argc, char *argv[]) {
//  // initialise the random number generator
//  srand(time(NULL));
//  
//  printf("Allocating memory for required data structures...");
//  int *arrow = malloc(N*sizeof(int));
//  if (arrow == NULL) {
//    printf("Unable to allocate requisite memory for program. Exiting..");
//    return -1;
//  }
//  
//  int *step = malloc(N*sizeof(int));
//  if (step == NULL) {
//    printf("Unable to allocate requisite memory for program. Exiting..");
//    return -1;
//  }
//
//  int *lasttime = malloc(N*sizeof(int));
//  if (lasttime == NULL) {
//    printf("Unable to allocate requisite memory for program. Exiting..");
//    return -1;
//  }
//
//  int iter;
//  
//  int **arrows = malloc(L*sizeof(int *));
//  
//  for (int i = 0; i < L; i++) {
//    arrows[i] = (int *)malloc(L*sizeof(int));
//  }
//
//  if (arrows == NULL) {
//    printf("Bad memory allocation. Unable to continue. Exiting..");
//    return -1;
//  }
//
//  struct vertex **vertices = malloc(L*sizeof(struct vertex *));
//
//  for (int i = 0; i < L; i++) {
//    vertices[i] = (struct vertex *)malloc(L*sizeof(struct vertex));
//  }
//
//  if (vertices == NULL) {
//    printf("Bad memory allocation. Unable to continue. Exiting..");
//    return -1;
//  }
//
//  generate_initial_lattice_state(arrows, vertices);
////
////  for (iter = 0; iter < 1000000; iter++) {
////    move(arrow, step, lasttime, iter);
////
////    if (iter % 5000) {
////      if (store_lattice_state(arrow) == -1) {
////        printf("Failed to write lattice state to file at step %d", iter);
////        continue;
////      }
////    }
////  }
////
//  free(arrows);
//  free(vertices);
//  return 0;
//}

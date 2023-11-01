/*
 * File: barnes_hut.c: Implements the Barnes Hut algorithm for n-body
 * simulation with galaxy-like initial conditions.
 */
#include "barnes_hut.h"

#define N 100000
#define NODE_SIZE sizeof(struct node_t)
#define CHILDREN_ARRAY_SIZE NODE_SIZE * 4

 //Some constants and global variables
const double L = 1, W = 1, dt = 1e-3, alpha = 0.25, V = 50, epsilon = 1e-1, grav = 0.04; //grav should be 100/N
double x[N], y[N], u[N], v[N], force_x[N], force_y[N];
double mass[N];
struct node_t* root;
char* file = "galaxy.in";

/*
  * Function to read a case
*/
int read_case(char* filename) {
  int i;
  FILE* arq = fopen(filename, "r");

  #pragma omp parallel for private(i)
  for (i = 0; i < N; i++) {
    fscanf(arq, "%lf %lf %lf %lf %lf", &mass[i], &x[i], &y[i], &u[i], &v[i]);
  }

  fclose(arq);
  return 0;
}

/*
 * Prints statistics: time, N, final velocity, final center of mass
 */
void print_statistics(double s, double e, float ut, float vt, float xc, float xy) {
  printf("Time elapsed in seconds: %f\n", e - s);
  printf("%.5f %.5f\n", ut, vt);
  printf("%.5f %.5f\n", xc, xy);
}

/*
 * Updates the positions of the particles of a time step.
 */
void time_step(void) {
  //Allocate memory for root
  root = malloc(NODE_SIZE);

  set_node(root);

  root->min_x = 0;
  root->max_x = 1;
  root->min_y = 0;
  root->max_y = 1;

  //Put particles in tree
  for (int i = 0; i < N; i++) {
    put_particle_in_tree(i, root);
  }

  //Calculate mass and center of mass
  calculate_mass(root);
  calculate_center_of_mass_x(root);
  calculate_center_of_mass_y(root);

  //Calculate forces
  update_forces();

  //Update velocities and positions
  double ax = 0.0, ay = 0.0;
  int i;
  #pragma omp parallel for private(i) reduction(+: ax, ay)
  for (i = 0; i < N; i++) {
    ax = force_x[i] / mass[i];
    ay = force_y[i] / mass[i];
    u[i] += ax * dt;
    v[i] += ay * dt;
    x[i] += u[i] * dt;
    y[i] += v[i] * dt;

    /* This of course doesn't make any sense physically,
     * but makes sure that the particles stay within the
     * bounds. Normally the particles won't leave the
     * area anyway.
     */
    bounce(&x[i], &y[i], &u[i], &v[i]);
  }

  //Free memory
  free_node(root);
  free(root);
}

/*
 * If a particle moves beyond any of the boundaries then bounce it back
 */
void bounce(double* x, double* y, double* u, double* v) {
  if (*x > 1.0f) {
    *x = 2 - *x;
    *u = -*u;
  }
  else if (*x < 0) {
    *x = -*x;
    *u = -*u;
  }

  if (*y > 1.0f) {
    *y = 2 - *y;
    *v = -*v;
  }
  else if (*y < 0) {
    *y = -*y;
    *v = -*v;
  }
}

/*
 * Puts a particle recursively in the Barnes Hut quad-tree.
 */
void put_particle_in_tree(int new_particle, struct node_t* node) {
  //If no particle is assigned to the node
  if (!node->has_particle) {
    node->particle = new_particle;
    node->has_particle = 1;
    return;
  }

  double half_current_x = (node->min_x + node->max_x) / 2;
  double half_current_y = (node->min_y + node->max_y) / 2;

  //If the node has no children
  if (!node->has_children) {
    //Allocate and initiate children
    node->children = malloc(CHILDREN_ARRAY_SIZE);
    struct node_t* children_0_add = &node->children[0];
    struct node_t* children_1_add = children_0_add + 1;
    struct node_t* children_2_add = children_0_add + 2;
    struct node_t* children_3_add = children_0_add + 3;

    set_node(children_0_add);
    set_node(children_1_add);
    set_node(children_2_add);
    set_node(children_3_add);

    //Set boundaries for the children
    children_0_add->min_x = node->min_x;
    children_0_add->max_x = half_current_x;
    children_0_add->min_y = node->min_y;
    children_0_add->max_y = half_current_y;

    children_1_add->min_x = half_current_x;
    children_1_add->max_x = node->max_x;
    children_1_add->min_y = node->min_y;
    children_1_add->max_y = half_current_y;

    children_2_add->min_x = node->min_x;
    children_2_add->max_x = half_current_x;
    children_2_add->min_y = half_current_y;
    children_2_add->max_y = node->max_y;

    children_3_add->min_x = half_current_x;
    children_3_add->max_x = node->max_x;
    children_3_add->min_y = half_current_y;
    children_3_add->max_y = node->max_y;

    //Put old particle into the appropriate child
    place_particle(node->particle, node, half_current_x, half_current_y);

    //It now has children
    node->has_children = 1;
  }

  //Add the new particle to the appropriate children
  //Put new particle into the appropriate child
  place_particle(new_particle, node, half_current_x, half_current_y);
}


/*
 * Puts a particle in the right child of a node with children.
 */
void place_particle(int particle, struct node_t* node, double half_current_x, double half_current_y) {
  double particle_x = x[particle];
  double particle_y = y[particle];

  if (particle_x <= half_current_x && particle_y <= half_current_y) {
    put_particle_in_tree(particle, &node->children[0]);
    return;
  }

  if (particle_x >= half_current_x && particle_y >= half_current_y) {
    put_particle_in_tree(particle, &node->children[3]);
    return;
  }

  if (particle_x < half_current_x && particle_y > half_current_y) {
    put_particle_in_tree(particle, &node->children[2]);
    return;
  }

  put_particle_in_tree(particle, &node->children[1]);
}

/*
  Sets initial values for a new node
*/
void set_node(struct node_t* node) {
  node->has_particle = 0;
  node->has_children = 0;
}

/*
  Frees memory for a node and its children recursively.
*/
void free_node(struct node_t* node) {
  if (node->has_children) {
    free_node(&node->children[0]);
    free_node(&node->children[1]);
    free_node(&node->children[2]);
    free_node(&node->children[3]);
    free(node->children);
  }
}

/*
  Calculates the total mass for the node. It recursively updates the mass
  of itself and all of its children.
*/
double calculate_mass(struct node_t* node) {
  if (!node->has_particle) {
    node->total_mass = 0;
  }
  else if (!node->has_children) {
    node->total_mass = mass[node->particle];
  }
  else {
    node->total_mass = 0;

    node->total_mass += calculate_mass(&node->children[0]);
    node->total_mass += calculate_mass(&node->children[1]);
    node->total_mass += calculate_mass(&node->children[2]);
    node->total_mass += calculate_mass(&node->children[3]);
  }

  return node->total_mass;
}

/*
  Calculates the x-position of the centre of mass for the
  node. It recursively updates the position of itself and
  all of its children.
*/
double calculate_center_of_mass_x(struct node_t* node) {
  if (!node->has_children) {
    node->c_x = x[node->particle];
  }
  else {
    node->c_x = 0;
    double m_tot = 0;

    if (node->children[0].has_particle) {
      node->c_x += node->children[0].total_mass * calculate_center_of_mass_x(&node->children[0]);
      m_tot += node->children[0].total_mass;
    }

    if (node->children[1].has_particle) {
      node->c_x += node->children[1].total_mass * calculate_center_of_mass_x(&node->children[1]);
      m_tot += node->children[1].total_mass;
    }

    if (node->children[2].has_particle) {
      node->c_x += node->children[2].total_mass * calculate_center_of_mass_x(&node->children[2]);
      m_tot += node->children[2].total_mass;
    }

    if (node->children[3].has_particle) {
      node->c_x += node->children[3].total_mass * calculate_center_of_mass_x(&node->children[3]);
      m_tot += node->children[3].total_mass;
    }

    node->c_x /= m_tot;
  }

  return node->c_x;
}

/*
  Calculates the y-position of the centre of mass for the
  node. It recursively updates the position of itself and
  all of its children.
*/
double calculate_center_of_mass_y(struct node_t* node) {
  if (!node->has_children) {
    node->c_y = y[node->particle];
  }
  else {
    node->c_y = 0;
    double m_tot = 0;

    if (node->children[0].has_particle) {
      node->c_y += node->children[0].total_mass * calculate_center_of_mass_y(&node->children[0]);
      m_tot += node->children[0].total_mass;
    }

    if (node->children[1].has_particle) {
      node->c_y += node->children[1].total_mass * calculate_center_of_mass_y(&node->children[1]);
      m_tot += node->children[1].total_mass;
    }

    if (node->children[2].has_particle) {
      node->c_y += node->children[2].total_mass * calculate_center_of_mass_y(&node->children[2]);
      m_tot += node->children[2].total_mass;
    }

    if (node->children[3].has_particle) {
      node->c_y += node->children[3].total_mass * calculate_center_of_mass_y(&node->children[3]);
      m_tot += node->children[3].total_mass;
    }

    node->c_y /= m_tot;
  }
  return node->c_y;
}

/*
  Calculates the forces in a time step of all particles in
  the simulation using the Barnes Hut quad tree.
*/
void update_forces() {
  int i;
  #pragma omp parallel for private(i) schedule(guided)
  for (i = 0; i < N; i++) {
    force_x[i] = 0;
    force_y[i] = 0;
    update_forces_help(i, root);
  }
}

/*
  Help function for calculating the forces recursively
  using the Barnes Hut quad tree.
*/
void update_forces_help(int particle, struct node_t* node) {
  //The node is a leaf node with a particle and not the particle itself
  if (node->has_children) {
    double aux = x[particle] - node->c_x;
    double temp = y[particle] - node->c_y;
    double r = sqrt(aux * aux + temp * temp);

    /* If the distance to the node's centre of mass is far enough, calculate the force,
     * otherwise traverse further down the tree
     */
    if ((node->max_x - node->min_x) / r < 0.5) {
      calculate_force(particle, node, r);
      return;
    }

    update_forces_help(particle, &node->children[0]);
    update_forces_help(particle, &node->children[1]);
    update_forces_help(particle, &node->children[2]);
    update_forces_help(particle, &node->children[3]);
  }
  //The node has children
  else if (node->particle != particle && node->has_particle) {
    //Calculate r 
    double aux = x[particle] - node->c_x;
    double temp = y[particle] - node->c_y;
    double r = sqrt(aux * aux + temp * temp);
    calculate_force(particle, node, r);
  }
}

/*
  Calculates and updates the force of a particle from a node.
*/
void calculate_force(int particle, struct node_t* node, double r) {
  double aux = r + epsilon;
  double temp = -grav * mass[particle] * node->total_mass / (aux * aux * aux);
  force_x[particle] += (x[particle] - node->c_x) * temp;
  force_y[particle] += (y[particle] - node->c_y) * temp;
}

int main(int argc, char* argv[]) {
  // The second argument sets the number of time steps
  int time_steps = atoi(argv[1]);

  if (read_case(file) == 1) {
    printf("Error: case instantiation failed.\n");
    return 1;
  }

  //Begin taking time
  double start = omp_get_wtime();

  //The main loop
  for (int i = 0; i < time_steps; i++) {
    time_step();
  }

  //Stop taking time
  double stop = omp_get_wtime();

  //Compute final statistics
  double vu = 0;
  double vv = 0;
  double sumx = 0;
  double sumy = 0;
  double total_mass = 0;
  int i;

  #pragma omp parallel for private(i) reduction(+:sumx,sumy,vu,vv,total_mass)
  for (i = 0; i < N; i++) {
    sumx += mass[i] * x[i];
    sumy += mass[i] * y[i];
    vu += u[i];
    vv += v[i];
    total_mass += mass[i];
  }

  double cx = sumx / total_mass;
  double cy = sumy / total_mass;

  print_statistics(start, stop, vu, vv, cx, cy);

  return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_STRLEN 256 // length of string for file name
#define COMMOUT    '#' // comment out
#define DELIMIT    '=' // delimiter between parameter tag and values
#define X_MIN      -5 // Minimum X for simulation field
#define X_MAX       5 // Maximum X for simulation field
#define Y_MIN      -5 // Minimum Y for simulation field
#define Y_MAX       5 // Maximum Y for simulation field
#define MAX_STEPS 10000 // Maximum time steps
#define BOUNDARY    1.0 // x position with maximum probability of piruette

/*
  unit of distance: mm
  unit of time: sec
*/

typedef
struct PARAMETERS {
  double velocity;   // velocity
  double delta_time; // step size
  double head_theta; // angle of perturbation
  double head_len;   // distance between tip and center of mass
  double alpha_w;    // parameter alpha of probability function of weathervane
  double beta_w;     // parameter beta of probability function of weathervane
  double alpha_p;    // parameter alpha of probability function of piruette
  double beta_p;     // parameter beta of probability function of piruette
  double angle_p_m;  // median of angles changed for piruette
  double angle_p_s;  // variance of angles changed for piruette
  double angle_w_m;  // median of angles changed for weathervane
  double angle_w_s;  // variance of angles changed for weathervane
  int    vd_index;   // index of side touching to the ground. +1 or -1.
  int    sampling;   // number of times used for judging threshold of weathervane
  double boundary;   // median of position of change for locomotion direction. unit: millimeter
  double r0;         // distance from the "boundary" when t=0.
} PARAMS;

typedef
struct PH_DIST {
  double x;   // x coordinate
  double pH;  // pH
} PHDIST;

/* Functions */
void read_parameters(char*); // readming values from parameter file
int piruette_possibility(double); // probability of piruette. 1:ON, 0:OFF.
int weathervane_possibility(double); // probability of weathervane depending on concetration
int initial_position_probability(double); // initial direction of movemnt
void simulation(char*); // perform nematode moving simulation
double acid_conc(double); // pH distribution function
double piruette_angle(double, int*); // direction of moving after piruette determined randomly
double weathervane_angle(double*, int *, double, int*); // direction of moving after weathervane
int read_pH_distribution(char *); // read pH distribution data from file

/* Parameters */
PARAMS parameters; // parameters read from configuration file
PHDIST *x_pH;      // array of position-pH

// command line arguments
// 1st: file name with parameters
int main(int argc, char *argv[])
{
  char params_filename[MAX_STRLEN]; // name of parameter file
  char result_filename[MAX_STRLEN]; // fil name for saving simulation results
  char result_filename_template[MAX_STRLEN]; // template string for file names saving simulation results
  char acid_filename[MAX_STRLEN];  // name of file saving concentration distribution and probability parameter

  int num;  // number of nematodes
  int i;
  int x_pH_data_lines; // number of lines of file of position-pH
  struct timespec seed; // seed of random number

  // for DEBUG
  double x, probr, probc, delta_x, pH;

  // initialization of random number
  clock_gettime(CLOCK_REALTIME, &seed); // obtain system clock
  srand(seed.tv_nsec); // initialize random number

  // checking number of command line arguments
  if(argc>4) {
    strcpy(params_filename, argv[1]);
    read_parameters(params_filename); // read from file
    strcpy(result_filename_template, argv[2]);
    strcpy(acid_filename, argv[3]);
    num = atoi(argv[4]);
  }
  else {
    printf("Usage\n");
    printf("./threshold_sim config_file saving_file pH-distribution_file number\n");
    printf("   config_file: save result of movement, coordinates and angles etc\n");
    printf("   pH-distribution_file: position-pH data file\n");
    printf("   number: number of nematodes\n");
    exit(1);
  }

  // reading data of position-pH
  x_pH_data_lines  =read_pH_distribution(acid_filename);

  // start simulation
  for(i=0; i<num; i++) {
    // generate file names of saving data for each nematode
    sprintf(result_filename, "%s.%05d.dat", result_filename_template, i);

    simulation(result_filename); // perform simulation
  }
  
  free(x_pH); // free memory for array of position-pH data
  
  return 0;
}

// read data from parameter file
void read_parameters(char *filename)
{
  char temp_line[MAX_STRLEN]; // string temporary used for reading data
  char temp_tag[MAX_STRLEN];  // string temporary used for saving tag part
  int delimit_position; // position of delimiter "="
  FILE *params_fp;    // file pointer for parameter file
  
  // read from parameter file
  params_fp = fopen(filename, "r");   // open the file

  // read from file
  while(fgets(temp_line, MAX_STRLEN, params_fp)) {
    if(!isspace(temp_line[0]) && temp_line[0] != COMMOUT) { // not read lines with 1st char of #, space and line break
      if(strncmp(temp_line, "velo", 4)==0) { // velocity
	parameters.velocity
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "delta_time", 8)==0) { // time step
	parameters.delta_time
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      
      else if(strncmp(temp_line, "head_t", 6)==0) { // perturbation of moving direction
	parameters.head_theta
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "head_l", 6)==0) { // distance between tip and center of mass
	parameters.head_len
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      
      else if(strncmp(temp_line, "alpha_w", 7)==0) { // parameter alpha of probability function of weathervane
	parameters.alpha_w
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "beta_w", 6)==0) { // parameter beta of probability function of weathervane
	parameters.beta_w
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "alpha_p", 7)==0) { // parameter alpha of probability function of piruette
	parameters.alpha_p
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "beta_p", 6)==0) { // parameter beta of probability function of piruette
	parameters.beta_p
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }

      else if(strncmp(temp_line, "angle_p_m", 9)==0) { // median of angles changed for piruette
	parameters.angle_p_m
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "angle_p_s", 9)==0) { // variance of angles changed for weathervane
	parameters.angle_p_s
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      
      else if(strncmp(temp_line, "angle_w_m", 9)==0) { // median of angles changed for weathervane
	parameters.angle_w_m
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      else if(strncmp(temp_line, "angle_w_s", 9)==0) { // variance of angles changed for weathervane
	parameters.angle_w_s
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }
      
      else if(strncmp(temp_line, "vd_index", 8)==0) { // index of side touching to the ground
	parameters.vd_index
	  =(int)strtol(index(temp_line, (int)DELIMIT)+1, NULL, 0); // convert string after "=" to int
      }

      else if(strncmp(temp_line, "sampling", 8)==0) { // number of times used for judging threshold of weathervane
	parameters.sampling
	  =(int)strtol(index(temp_line, (int)DELIMIT)+1, NULL, 0); // convert string after "=" to int
      }

      else if(strncmp(temp_line, "boundary", 8)==0) { // median of position of change for locomotion direction
	parameters.boundary
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }

      else if(strncmp(temp_line, "r0", 2)==0) { // distance from the "boundary" when t=0.
	parameters.r0
	  =strtod(index(temp_line, (int)DELIMIT)+1, NULL); // convert string after "=" to double
      }

      else {  // undefined tag
	delimit_position = (int)(strlen(temp_line)-strlen( index(temp_line, (int)DELIMIT) ) );
	temp_line[delimit_position]='\0';
	strncpy(temp_tag, temp_line, delimit_position+1);
	printf("Error in configuration file\n");
	printf("undefined tag => *** %s ***\n", temp_tag);
	
	fclose(params_fp); // close file
	
	exit(2); // force exiting
      }
    }
  }
  
  fclose(params_fp); // close file
}

// perform simulation
void simulation(char *outfilename)
{
  int i, j;
  double x00, y00; // initial position of center of mass (c.o.m)
  double x0, y0;   // coordinate of c.o.m
  double x1, y1;   // coordinate of tip
  double init_bound_xmin; // minimun x for initial position
  double vx, vy;   // vector of moving direction of c.o.m
  double delta_angle; // direction change of c.o.m, angle between x axis
  double zeta;     // direction of moving of c.o.m, angle between x axis
  int p_piruette;   // piruette or NOT: 0:NO, 1:YES
  int p_weathervane;// weathervane or NOT: 0:NO, 1:YES
  int steps;       // total step number of simulation
  double conc[2];  // concentration at the tip
                   // conc[0]: conc. at previous step. conc[1]: conc at current step.
  double delta_angle_pir; // angle change by piruette
  double delta_angle_wv;  // angle change by weathervane
  double cos_steps[2]; // status of tip position
                       // [0]: position at previous step. [1]: at current step.
  int wv_flag =0; // flag for weathervane or not at the next step
  FILE *result_fp;

  // simulation
  steps =-1; // count time step during simulation
  
  // initial x
  init_bound_xmin = ((-parameters.r0 +parameters.boundary > X_MIN ) ?
		     -parameters.r0 +parameters.boundary : X_MIN);

  // initial poistion (tail)
  x0 = ((double)rand()/RAND_MAX)*(parameters.boundary - init_bound_xmin) + init_bound_xmin;
  y0 = ((double)rand()/RAND_MAX)*(Y_MAX - Y_MIN) + Y_MIN;

  // generation of initial moving direction
  do {
    zeta = ((double)rand()/RAND_MAX-0.5)*M_PI;
  } while( initial_position_probability(zeta)!=1 );  // adopt x0,y0 stochastically according to sin(zeta)
  
  // initial velocity
  vx = parameters.velocity *cos(zeta);
  vy = parameters.velocity *sin(zeta);

  // initial position of tip. index of head direction.
  cos_steps[0] = cos(M_PI/2*steps); // the index changes in order of +1, 0, -1, 0.
  x1 = x0 + parameters.head_len *cos(zeta +parameters.head_theta/180.0 *M_PI*cos_steps[0]);
  y1 = y0 + parameters.head_len *sin(zeta +parameters.head_theta/180.0 *M_PI*cos_steps[0]);

  conc[1] = acid_conc(x1); // concentration at the current position
  
  // file for saving results
  result_fp = fopen(outfilename, "w");
  
  // initialize angle change by piruette and weathervane
  delta_angle_pir = 0.0;
  delta_angle_wv  = 0.0;

  // initialize probalitities of piruette and weathervane. dummy.
  p_piruette = 0;
  p_weathervane = 0;

  // save information
  fprintf(result_fp, "%d %lf %lf %lf %lf %lf %lf %d %lf %d %lf\n",
	  steps, conc[1], x0, y0, x1, y1, zeta,
	  p_piruette, delta_angle_pir, p_weathervane, delta_angle_wv);
  
  steps++;

  while(steps <MAX_STEPS &&
	x0>X_MIN && x0<X_MAX &&
	y0>Y_MIN && y0<Y_MAX) {

    cos_steps[1] = cos(M_PI/2*steps);

    // measuring concentration
    // concentration at previous time steop
    conc[1] = acid_conc(x1);
    
    // change of moving direction by piruette
    delta_angle_pir =  parameters.vd_index * piruette_angle(conc[1], &p_piruette);
    
    // change of moving direction by weathervane
    if(p_piruette != 1) { // weathrevane occurs when piruette does not occur
      delta_angle_wv = weathervane_angle(conc, &wv_flag, cos_steps[1]-cos_steps[0], &p_weathervane);
    }
    else {  // weathrevane does not occur when piruette occur
      delta_angle_wv = 0.0;
      wv_flag =0;  // reset status of weathervane
      p_weathervane = 0; // make probability of weathervane zero
    }

    // change of moving direction (angle between x axis)
    delta_angle = delta_angle_pir + delta_angle_wv;
    
    // new moving direction
    zeta += delta_angle;
      
    // renew current concentration
    conc[0] = conc[1];
    // end of measuring concentraion

    // velocity vector of c.o.m
    vx = parameters.velocity * cos(zeta);
    vy = parameters.velocity * sin(zeta);

    // position of c.o.m
    x0 = x0 + vx*parameters.delta_time;
    y0 = y0 + vy*parameters.delta_time;

    // position of tip
    x1 = x0 + parameters.head_len *cos(zeta +parameters.head_theta/180.0*M_PI *cos_steps[1]);
    y1 = y0 + parameters.head_len *sin(zeta +parameters.head_theta/180.0*M_PI *cos_steps[1]);

    // save information
    fprintf(result_fp, "%d %lf %lf %lf %lf %lf %lf %d %lf %d %lf\n",
	    steps, conc[1], x0, y0, x1, y1, zeta,
	    p_piruette, delta_angle_pir, p_weathervane, delta_angle_wv);
    
    // reset change of angle by piruette and weathervane
    delta_angle_pir = 0.0;
    delta_angle_wv  = 0.0;

    cos_steps[0] = cos_steps[1];    // save current index of tip direction
    
    steps++; // count step
  }
  
  fclose(result_fp); // close file
}

// random change (absolute value) of moving direction after piruette
double piruette_angle(double conc, int *p_piruette)
{
  double x, y;     // Gaussian random number
  double r1, r2;   // uniform random number
  int sign;        // index for clockwise or counterclockwise
  double turning_angle; // angle of change for moving direction. absolute value.

  *p_piruette = piruette_possibility(conc); // probability of piruette

  if ( !(*p_piruette) ) {// when p_piruette=0, piruette does not occur.
    return 0.0;
  }
  
  // generation of Gaussian random number by Box-Muller method
  r1 = (double)rand()/RAND_MAX;
  r2 = (double)rand()/RAND_MAX;
  x = parameters.angle_p_s *sqrt(-2*log(r1)) *cos(2*M_PI*r2)
    +parameters.angle_p_m;

  // absolute value of angle
  turning_angle = x/180.0*M_PI; // convert degree to radian
  
  return turning_angle;
}

// change of moving direction after weathervane
double weathervane_angle(double *conc, int *flag, double cosindex, int *p_weathervane)
{
  double x, y;     // Gaussian random number
  double r1, r2;   // uniform random number
  int sign;        // index of clockwise or counterclockwise
  double turning_angle; // change of moving direction

  if(*flag) { // probability of weathervane is higher than threshold at the previous time step
    // generate Gaussian random number by Box-Muller method
    r1 = (double)rand()/RAND_MAX;
    r2 = (double)rand()/RAND_MAX;
    x = parameters.angle_w_s *sqrt(-2*log(r1)) *cos(2*M_PI*r2)
      +parameters.angle_w_m;

    // angle
    if(cosindex<0) { // current tip position is right side of the previous tip position
      if(conc[0] > conc[1] ) { // acid conc at the previous time step is higher than that at current step
	turning_angle = x/180.0*M_PI; // convert degree to radian
      }
      else {                  // acid conc at the previous time step is lower than that at current step
	turning_angle = -x/180.0*M_PI; // convert degree to radian
      }
    }
    else {  // current tip position is left side of the previous tip position
      if(conc[0] > conc[1] ) { // acid conc at the previous time step is higher than that at current step
	turning_angle = -x/180.0*M_PI; // convert degree to radian
      }
      else {                  // acid conc at the previous time step is lower than that at current step
	turning_angle = x/180.0*M_PI; // convert degree to radian
      }
    }
    
    // resetting flag
    *flag = 0;
  
    return turning_angle;
  }
  else { // probability of weathervane is below threshold
    *p_weathervane = weathervane_possibility(conc[1]); // probability of weathervane

    if( (*p_weathervane) ) { // when p_weathervane=1
      *flag = 1; // flag indicating "weathervane is occuring"
    }
    return 0.0;  // at this line, weathervane has not occured
  }
}

// pH at position x
double acid_conc(double x)
{
  int i;
  double conc; // concentration. used as return value.

  // x_pH[i-1].x < x <= x_pH[i].x
  if(x>=X_MAX) x=X_MAX;
  if(x<=X_MIN) x=X_MIN;
  
  i=0;
  while(x_pH[i].x<x) {
    i++;
  }

  // obtain pH by interpolation
  conc
    = (x_pH[i].pH-x_pH[i-1].pH)/(x_pH[i].x-x_pH[i-1].x)*(x-x_pH[i-1].x)
    +x_pH[i-1].pH;
  
  return conc;
}

// probability of piruette depending on concentration
int piruette_possibility(double conc)
{
  double probability;
  double threshold;
  int piruette_idx;     // piruette occur: 1, not: 0

  threshold = (double)rand()/RAND_MAX; // uniform random number between 0 and 1
  probability = 1.0/( 1+exp((conc-parameters.beta_p)/parameters.alpha_p) );

  if(threshold < probability) {
    piruette_idx = 1;
  }
  else {
    piruette_idx = 0;
  }

  return piruette_idx;
}

// probability of weathervane depending on concentration
int weathervane_possibility(double conc)
{
  double probability;
  double threshold;
  int weathervane_idx;     // weathervane occur: 1, not: 0

  threshold = (double)rand()/RAND_MAX; // random number
  probability = 1.0/( 1+exp((conc-parameters.beta_w)/parameters.alpha_w) );

  if(threshold < probability) {
    weathervane_idx = 1;
  }
  else {
    weathervane_idx = 0;
  }

  return weathervane_idx;
}

// probability of adopting initial position depending of initial moving direction
int initial_position_probability(double zeta)
{
  double probability;
  double threshold;
  int position_idx;     // adopt: 1, not: 0

  threshold = (double)rand()/RAND_MAX; // random number
  probability = cos(zeta);

  if(threshold < probability) {
    position_idx = 1;
  }
  else {
    position_idx = 0;
  }

  return position_idx;
}

// read pH distribution data
int read_pH_distribution(char *acid_filename)
{
  int lines; // counter for line number
  int i;
  char temp_line[MAX_STRLEN]; // temporal string for reading file
  char *tp;  // string for dividing from longer string
  char *saveptr; // pointer for strtok_r
  FILE *distr_fp; // file pointer

  // open file
  distr_fp = fopen(acid_filename, "r");

  // obtain line number in the file
  lines = 0;
  while(fgets(temp_line, MAX_STRLEN, distr_fp)) {
    lines++;
  }

  // rewind to top of the file
  rewind(distr_fp);

  // allocate array of structure with the number of lines in the file
  // array of structure x_pH is global
  x_pH=(PHDIST *)malloc(sizeof(PHDIST) *lines);

  // read position and pH from file and store to the array of structure
  i=0;
  while(fgets(temp_line, MAX_STRLEN, distr_fp)) {
    tp =strtok_r(temp_line, " ", &saveptr);
    x_pH[i].x = strtod(tp, NULL); // position at the first line
    
    while( tp != NULL ){ // read until empty line
      tp =strtok_r(NULL, " ", &saveptr);
      if(tp != NULL ){
	x_pH[i].pH = strtod(tp, NULL); // pH stored in 2nd column
      }
    }
    i++;
  }

  fclose(distr_fp); // close file

  return lines; // number of lines is return value
}

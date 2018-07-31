#include "ranlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "mpi.h"
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std ;

// dim : number of species
// R :   number of reactions
extern int dim, R ;

// num_basis : number of basis functions
extern int num_basis ;

// N_traj : total number of trajectories

/* local_N_traj : number of trajectories that are distributed
 *  		  on a local process */
extern int mpi_rank , mpi_size , N_traj, local_N_traj, local_traj_start_idx ;

/* know_reactions_flag:   equals 1 if reaction types are known, 
 			  otherwise equals 0 */

// xx_basis_flag : 	  equals 1 if x*x is used as basis, 
//       	   	  equals 0 if x*(x-1) is used
extern int know_reactions_flag, xx_basis_flag ;

// num_state_in_traj : number of states in each trajectory 
extern vector<int> num_state_in_traj ;

// total number of reactions in all N_traj trajectories
extern int tot_num_reaction_in_traj ;

// numerical scheme used to solve the (sparse) optimization problem
extern int iter_scheme ;

// T :   length of trajectories when they are simulated using SSA method
extern double T ;

// total length of time for all N_traj trajectories
extern double total_T ;

// T_traj_vec : vector containing length of time for each trajectory 
extern vector<double> T_traj_vec; 

// sparse intensity constant, the sparse intensity for each 
// parameter (with indice i, j) is regular_lambda * omega_weights[i][j]
extern double regular_lambda ;

// delta :      constant \epsilon used in the function G_\epsilon
//
// stop_eps :   iteration schemes stop when the difference of two steps is not
// 	        larger than this value
//
// descent_dt : the step-size used in the iterative scheme
//
// g_cut :      function G(x) returns x when x/delta >= g_cut
//
// eps   : 	parameter used in the epsL1_norm (an approximation of l^1 norm)
extern double stop_eps , eps , descent_dt, delta, g_cut ;

// tot_step : 		total iteration steps 
// output_interval : 	determine how often to print information during iteration 
extern int tot_step, output_interval ;

// maximal order of polynomial functions used as basis functions, in the
// current implementation we use poly_order = 1 or 2.
extern int poly_order ;

// init_state :     initial state of the trajectories in SSA method
// reactant_num :   number of reactants in each reactions
extern vector<int> init_state , reactant_num ; 

// Mi_in_traj : number of occurrences of reactions in channel i within the trajectory data
extern vector<vector<int> > Mi_in_traj ;

// Mi_in_all_traj : number of occurrences of reactions in channel i within all N_traj trajectories 
extern vector<int> Mi_in_all_traj ;

// in, out, and change vectors of each reactions
extern vector<vector<int> > vvec_in, vvec_out, vvec , v_neg_idx_vec ;

// true parameters, i.e., the constants in the propensity functions
extern vector<double> kappa_vec ;

// total number of unknown parameters need to be estimated
extern int total_unknown_omega_parameters ;

// unknown parameters and their sparsity weights  
extern vector<vector<double> > omega_vec , omega_weights ;

// basis_vec : 			basis functions used for reconstruction
// basis_index_per_channel :	indices of basis functions used for each channel 
extern vector<vector<int> > basis_vec, basis_index_per_channel ;

// waiting_time_vec :   vector of waiting times for each state in each trajectory
// t_vec : 		vector of the current time for each state in each trajectory
extern vector<vector<double> > waiting_time_vec, t_vec ;

// number of channels obtained from trajectories data
extern int channel_num ;

// channel_list : 		change vectors of each channel 
// reactions_in_channel : 	indices of channel for each reactions 
extern vector<vector<int> > channel_list, reactions_in_channel ;

// index of channels for each reaction within trajectory data
extern vector<vector<int> > channel_idx_in_traj ;

// vector of trajectories 
extern vector<vector<vector<int> > > traj_vec ;

// file pointer used to print log information 
extern FILE * log_file ;

int dir_check( char dir_name[] ) ;
int init(char * log_file_name) ;
double ssa(double & t_now, vector<int> & c_state, vector<int> & next_state) ;

double val_basis_funct(int basis_idx, vector<int> &state) ;
double shrinkage(double x, double lambda) ;
double grad_log_likelihood(vector<vector<double> > & coeff_vec, vector<vector<double> > & grad_coeff) ;
double log_likelihood(vector<vector<double> > & coeff_vec) ;
double val_ai(int channel_idx, vector<int> &state, vector<vector<double> > & coeff_vec ) ;
double val_basis_funct(int basis_idx, vector<int> &state) ;
void print_grad( vector<vector<double> > & coeff_vec ) ;
double l1_norm( vector<vector<double> > & coeff_vec, vector<vector<double> > & weights ) ;
double epsL1_norm( vector<vector<double> > & coeff_vec, vector<vector<double> > & weights ) ;

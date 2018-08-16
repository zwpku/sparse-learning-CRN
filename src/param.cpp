#include "sparse_learning.h"

int n , R , N_traj ;
int mpi_rank, mpi_size, local_N_traj, local_traj_start_idx ;

double T , total_T ;

vector<double> T_traj_vec ; 
double regular_lambda ;

double stop_eps , eps, delta, g_cut ;
int know_reactions_flag, xx_basis_flag ;
vector<int> num_state_in_traj ;

int tot_num_reaction_in_traj ;

int tot_step, output_interval ;

vector<vector<int> > Mi_in_traj ; 
vector<int> init_state , reactant_num , Mi_in_all_traj ;

vector<vector<int> > vvec_in, vvec_out, vvec, v_neg_idx_vec ;
vector<double> kappa_vec ;

int poly_order, num_basis, epsL1_flag, solver_id ;

double grad_dt ;

vector<vector<int> > basis_vec, basis_index_per_channel ;

int total_unknown_omega_parameters ;

vector<vector<double> > omega_vec, omega_weights ;

int channel_num ;

vector<vector<vector<int> > > traj_vec ;
vector<vector<int> > channel_list, reactions_in_channel ;
vector<vector<int> > channel_idx_in_traj ;

vector<vector<double> > waiting_time_vec, t_vec ;

FILE * log_file ;

#include "ssa.h"

// smooth rectifier G(x) = delta * ln(1+exp(x/delta)) 
double G(double x)
{
  if (know_reactions_flag == 1) return x ;
  else {
    if (x >= g_cut * delta) return x;
    if (x <= -g_cut * delta) return 0.0 ;
    return delta * log( 1 + exp( x / delta ) ) ;
  }
}

// derivative of G
double dG(double x)
{
  if (know_reactions_flag == 1) return 1 ;
  else 
  {
    if (x >= g_cut * delta) return 1.0 ;
    if (x <= -g_cut * delta) return 0.0 ;
    return 1.0 / (1 + exp(-x / delta)) ;
  }
}

// ln G(x)
double logG(double x)
{
  if (know_reactions_flag == 1) 
  {
    if (x <= 0) 
    {
      printf("Error: Can not compute ln(%.4f)!\n", x) ;
      exit(1) ;
    }
    return log(x) ;
  }
  else 
  {
    double tmp ;
    tmp = x / delta ;
    if (tmp > g_cut) 
      return log(x) + exp(-tmp) / tmp ;
    if (tmp < -g_cut)
      return x / delta + log(delta) ;
    return log(delta) + log( log(1 + exp(tmp)) ) ;
  }
}

// derivative of function ln G(x)
double d_logG(double x)
{
  if (know_reactions_flag == 1) 
  {
    if (x <= 0) 
    {
      printf("Error: x=%.4f <= 0! \n", x) ;
      exit(1) ;
    }
    return 1.0 / x;
  }
  else 
  {
    double tmp ;
    tmp = x / delta ;
    if (tmp > g_cut) return 1.0 / x ;
    if (tmp < -g_cut) return 1.0 / ( (1.0 + exp(tmp)) * delta ) ;
    return 1.0 / (delta * (1 + exp(-tmp)) * log(1 + exp(tmp))) ;
  }
}

// compute weighted l^1 norm of vector 
double l1_norm( vector<vector<double> > & coeff_vec, vector<vector<double> > & weights )
{
  double s;
  s = 0.0 ;
  for (int i =0 ; i < channel_num; i ++)
    for (int j = 0 ; j < coeff_vec[i].size() ; j ++)
    {
      s += fabs(coeff_vec[i][j]) * weights[i][j] ;
    }
  return s ;
}

// the absolute value |x| is approximated by (|x|^2+\epsilon)^{1/2}
double epsL1_norm( vector<vector<double> > & coeff_vec, vector<vector<double> > & weights )
{
  double s ;
  s = 0.0 ;
  for (int i =0 ; i < channel_num ; i ++)
    for (int j = 0 ; j < coeff_vec[i].size() ; j ++)
    {
      s += sqrt( coeff_vec[i][j] * coeff_vec[i][j] + eps ) * weights[i][j] ;
    }
  return s ;
}

/* 
 test whether the reaction with index r_idx is possible to happen, i.e., if there are enough
 reactants , then return 1, otherwise return 0
*/ 
int try_reaction( vector<int> & in_vec, int r_idx )
{
  for (int i = 0 ; i < dim; i ++)
  {
    if ( in_vec[i] < vvec_in[r_idx][i] ) return 0 ;
  }

  return 1 ;
}

// state : state vector
// r_idx : reaction id 
double rate_of_state( vector<int> & state, int r_idx )
{
  double tmp ;

  tmp = 1.0 ;
  for ( int j = 0 ; j < dim ; j ++ )
    for ( int k = 0 ; k < vvec_in[r_idx][j] ; k ++ )
      tmp *= ( state[j] - k ) * 1.0 ;

  return tmp ;
}

// compute the propensity function a_0 of state, and return the list of
// possible reactions 
double ssa_a( vector<int> & state , vector<int> & possible_reactions , vector<double> & a_vec ) 
{
  double a0 , tmp ;
  a0 = 0 ;

  possible_reactions.clear() ; 
  a_vec.clear() ;

  for ( int i = 0 ; i < R ; i ++ )
    if ( try_reaction(state, i) )
    {
      possible_reactions.push_back(i) ;
      tmp = kappa_vec[i] * rate_of_state(state, i) ;
      a_vec.push_back( tmp ) ;
      a0 += tmp ;
    }
  return a0 ;
}

// update the state of the system using SSA method
double ssa(double & t_now, vector<int> & c_state, vector<int> & next_state) 
{
  double tmp_r1, tmp_r2 , tmp , tau , a0 ;
  vector<int> possible_reactions ;
  vector<double> a_vec ;

  a0 = ssa_a( c_state, possible_reactions, a_vec ) ;

  tmp_r1 = ranf() ; 
  tmp_r2 = ranf() * a0 ; 

  // compute waiting time
  if (a0 > 0)
    tau = - 1.0 / a0 * log(tmp_r1) ;
  else  // no reaction can occur
    tau = T- t_now + 1 ;

  if (t_now + tau <= T) 
  {
    // select the reaction channel 
    tmp = 0 ;
    for (int j = 0 ; j < possible_reactions.size() ; j ++)
    {
      tmp += a_vec[j] ;
      if (tmp > tmp_r2)
      {
	// update the state of system by adding the change vector 
	for (int k = 0 ; k < dim ; k ++)
	  next_state[k] = c_state[k] + vvec[possible_reactions[j]][k] ;

	break ;
      }
    }
  }

  return tau ;
}

// compute the value of basis function
double val_basis_funct(int basis_idx, vector<int> &state)
{
  double s ;
  int idx1, idx2 ;

  if (basis_vec[basis_idx].size() == 1) // linear basis 
  {
    idx1 = basis_vec[basis_idx][0] ;
    s = state[idx1] ;
  }
  else  // quadratic basis 
  {
    idx1 = basis_vec[basis_idx][0] ;
    idx2 = basis_vec[basis_idx][1] ;

    if ( (idx1 == idx2) && (xx_basis_flag == 0) ) // x(x-1)
      s = state[idx1] * (state[idx1] - 1) ;
    else // x*y or x*x
      s = state[idx1] * state[idx2] ;
  }

  return s ;
}

// compute propensity ai of a reaction channel based on basis functions
double val_ai(int channel_idx, vector<int> &state, vector<vector<double> > & coeff_vec )
{
  int idx ;
  double s ; 
  s = 0 ; 
  for ( int i = 0 ; i < coeff_vec[channel_idx].size() ; i ++ )
  {
    idx = basis_index_per_channel[channel_idx][i] ;
    s += val_basis_funct(idx, state) * coeff_vec[channel_idx][i] ;
  }
  return s ;
}

// print some debug information
void dump_info(int idx, double tmp_ai, vector<int> & c_state, vector<vector<double> > & coeff_vec )
{
  int ii ;
  printf("\t State : ");
  for (int ii = 0 ; ii < dim ; ii ++)
    printf("%d ", c_state[ii]); 
  printf("\n\t Channel idx: %d\n", idx) ;
  printf("\t Omega coefficients: ");
  for (int ii = 0 ; ii < coeff_vec[idx].size(); ii ++)
    printf("%.4lf\t", coeff_vec[idx][ii]); 
  printf("\n\t ai = %.4lf\n", tmp_ai);
}

// compute log-likelihood function
double log_likelihood(vector<vector<double> > & coeff_vec)
{
  double s, tmp_ai, s1, ai, local_s ;
  int idx ;

  local_s = 0 ;
  // process trajectories that are distributed on the local processor
  for (int traj_idx = 0; traj_idx < local_N_traj; traj_idx ++)
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++) // loop for each jump (or reaction)
  {
    // subtract by 1, because index starts from one   
    idx = channel_idx_in_traj[traj_idx][i] - 1 ;
    tmp_ai = val_ai(idx, traj_vec[traj_idx][i], coeff_vec) ;

    // since the channel occurs in trajectory, the propensity ai should be
    // positive
    if (tmp_ai == 0)
    {
      printf( "Error: ai=0! i=%d, idx=%d", i, idx ) ;
      fprintf( log_file, "Error: ai=0! i=%d, idx=%d", i, idx ) ;

      dump_info(idx, tmp_ai, traj_vec[traj_idx][i], coeff_vec ) ;
      exit(1) ;
    }

    // compute the first part of the log-likelihood function
    local_s += logG(tmp_ai) ;
  }

  // process trajectories that are distributed on the local processor
  for (int traj_idx = 0; traj_idx < local_N_traj; traj_idx ++)
  // loop for each state, including the last one
    for (int i = 0; i < num_state_in_traj[traj_idx]; i ++)
    {
      s1 = 0.0 ;
      for (int j = 0 ; j < channel_num; j ++)
      {
	tmp_ai = val_ai( j, traj_vec[traj_idx][i], coeff_vec ) ;
	s1 += G(tmp_ai) ;
      }
    // compute the second part of the log-likelihood function
      local_s -= waiting_time_vec[traj_idx][i] * s1 ;
    }

  // sum up all processors
  MPI_Allreduce(&local_s, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;

  return s ;
}

double grad_log_likelihood(vector<vector<double> > & coeff_vec, vector<vector<double> > & grad_coeff)
{
  double s, tmp_ai, s1, tmp , local_s ;
  int idx , basis_idx ;

  // compute each component, indexed by i0 (channel) and j0 (basis function)
  for (int i0 = 0 ; i0 < channel_num ; i0 ++)
    for (int j0 = 0 ; j0 < coeff_vec[i0].size() ; j0 ++)
    {
      // index of basis function
      basis_idx = basis_index_per_channel[i0][j0] ;
      // compute the gradient of the first (log) part 
      local_s = 0.0; 

      // loop for each trajectory on the local processor
      for (int traj_idx = 0; traj_idx < local_N_traj; traj_idx ++)
      // loop for each jump (reaction)
	for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++)
	{
	// subtract by 1, because index starts from one   
	  idx = channel_idx_in_traj[traj_idx][i] - 1 ;
	  // only when the jump is in channel i0 
	  if (idx != i0) continue ;
	  // derivative by chain rule
	  tmp = val_basis_funct( basis_idx, traj_vec[traj_idx][i] ) ;
	  tmp_ai = val_ai(idx, traj_vec[traj_idx][i], coeff_vec ) ;
	  local_s -= tmp * d_logG(tmp_ai) ;
	}

      // gradient of the second part
      for (int traj_idx = 0; traj_idx < local_N_traj; traj_idx ++)
	for (int i = 0; i < num_state_in_traj[traj_idx]; i ++)
	{
	  tmp = val_basis_funct( basis_idx, traj_vec[traj_idx][i] ) ;
	  tmp_ai = val_ai(i0, traj_vec[traj_idx][i], coeff_vec ) ;
	  local_s += tmp * dG(tmp_ai) * waiting_time_vec[traj_idx][i] ;
	}

      // sum up all trajectories among different processors
      MPI_Allreduce(&local_s, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
      
      // normalize by dividing the total time 
      grad_coeff[i0][j0] = s / total_T ;
    }
}

// shrinkage operator: max(|x|-lambda, 0) * sgn(x)
double shrinkage(double x, double lambda)
{
  if (x > lambda) 
    return x - lambda ;
  if (x < -lambda)
    return x + lambda ;
  return 0.0;
}

void print_grad( vector<vector<double> > & coeff_vec )
{
  for (int i =0 ; i < channel_num; i ++)
  {
    printf("Channel %d (", i) ;
    for (int j = 0 ; j < channel_list[i].size() ; j ++)
      cout << channel_list[i][j] << ' ';
    cout << ") :\t";
    for (int j = 0 ; j < coeff_vec[i].size() ; j ++)
    {
      cout << coeff_vec[i][j] << '\t' ;
    }
    cout << endl ;
  }
}

// check whether a directory exists or not 
int dir_check( char dir_name[] )
{
  struct stat sd ;

  if ( (stat(dir_name, &sd) != 0) || (S_ISDIR(sd.st_mode) != 1) )
    {
      printf( "\nNotice: directory doesn't exist : %s\n", dir_name ) ;
      // try to create a new one
      if ( mkdir(dir_name, 0775) == 0 ) 
	printf( "\nNotice: directory %s is created. \n\n", dir_name ) ;
      else return -1 ;
    }
  return 0 ;
}


#include "sparse_learning.h"

/*
 *
 * Smooth rectifier G_\delta(x) = \delta * ln(1+exp(x/\delta)) 
 *
 * G_\delta(x) is a smooth approximation of the activation function \max(x,0),
 * since we know
 *
 * 	\lim_{\delta \rightarrow 0} G_\delta(x) = \max(x,0)
 *
 * g_cut:   	parameter used to avoid numerical error when x is too large or
 * 		too small.
 */
double G(double x)
{
  if (know_reactions_flag == 1) return x ;
  else {
    if (x >= g_cut * delta) return x;
    if (x <= -g_cut * delta) return 0.0 ;
    return delta * log( 1 + exp( x / delta ) ) ;
  }
}

/* 
 * compute 1st derivative of G_\delta(x)
 *
 */
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

/* 
 * Compute ln G(x)
 */
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

/* 
 * compute the 1st derivative of the function ln G(x)
 */
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

/*
 * check whether x equals zero
 */
int is_zero(double x)
{
  if (fabs(x) < 1e-15) return 1 ;
  else return 0 ;
}

/* 
 * Compute the weighted l^1 norm of vector (corresponding to certain reaction
 * channel
 *
 * input :
 *   i : 		index of the channel
 *   coeff_vec :	coefficient vector
 *   weights :		the weights for each component
 *
 * return :
 *   s :		the wegithed l^1 norm
 */
double l1_norm_partial( int i, vector<vector<double> > & coeff_vec, vector<vector<double> > & weights )
{
  double s;
  s = 0.0 ;
    for (int j = 0 ; j < coeff_vec[i].size() ; j ++)
    {
      s += fabs(coeff_vec[i][j]) * weights[i][j] ;
    }
  return s ;
}

/* 
 * approximation of the weighted l^1 norm, where
 * the absolute value |x| is approximated by (|x|^2+\epsilon)^{1/2}
 */
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
 * Test whether a reaction is possible to happen
 *
 * input :
 *   in_vec :	current state of the system
 *   r_idx :	index of the reaction
 *
 * return :  
 *   1, 	if the reaction can happen, i.e., there are enough reactants 
 *   0, 	otherwise
 *
*/ 
int try_reaction( vector<int> & in_vec, int r_idx )
{
  for (int i = 0 ; i < n ; i ++)
  {
    if ( in_vec[i] < vvec_in[r_idx][i] ) return 0 ;
  }

  return 1 ;
}

/*
 * propensity of a given state 
 *
 * input:
 * 	state: 	current state of the system
 * 	r_idx: 	index of reaction 
 *
 * output: 
 *   	Essentially, the function compute either 
 *    		x, x*y, or x(x-1),
 *   	depending on the state and the type of the reaction.
 *   	See the document in the file prepare.cpp for details.
 *
 */
double rate_of_state( vector<int> & state, int r_idx )
{
  double tmp ;

  tmp = 1.0 ;
  for ( int j = 0 ; j < n ; j ++ )
    for ( int k = 0 ; k < vvec_in[r_idx][j] ; k ++ )
      tmp *= ( state[j] - k ) * 1.0 ;

  return tmp ;
}

/* 
 *
 * compute the (total) propensity function a_0 of state, and return the list of
 * possible reactions 
 *
 * input:
 *   state: 			current state of the system
 *   possible_reactions: 	list of indices of reactions that can occur at
 *   				current state
 *   a_vec:			propensity functions of each reaction that can
 *   				occur
 *
 * return: 
 *   a0: 			the (total) propensity function (equals to the
 *   				sum of numbers in the vector a_vec)
 *
 */

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

/*
 *
 * update the state of the system using SSA method
 *
 * input:
 *   t_now : 		current time
 *   c_state : 		current state
 *   next_state : 	new state after the next reaction occurs
 *
 * return : 
 *   tau : 		the waiting time before the next reaction occurs
 *
 */
double ssa(double & t_now, vector<int> & c_state, vector<int> & next_state) 
{
  double tmp_r1, tmp_r2 , tmp , tau , a0 ;
  vector<int> possible_reactions ;
  vector<double> a_vec ;

  // compute the (total) propensity 
  a0 = ssa_a( c_state, possible_reactions, a_vec ) ;

  // generate two random variables
  tmp_r1 = ranf() ; 
  tmp_r2 = ranf() * a0 ; 

  // compute waiting time
  if (a0 > 0) // exponential distribution 
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
	for (int k = 0 ; k < n ; k ++)
	  next_state[k] = c_state[k] + vvec[possible_reactions[j]][k] ;

	break ;
      }
    }
  }

  return tau ;
}

/* 
 *
 * compute the value of basis function
 *
 * input:
 *   basis_idx :	index of the basis function
 *   state :		current state of the system
 *
 * output: 
 *   s : 		value of the basis function
 *
 */
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

    /* 
     *
     * When xx_basis_flag=0, the basis function is x*(x-1)
     *
     * When xx_basis_flag=1, the basis function is x*x
     *
     */
    if ( (idx1 == idx2) && (xx_basis_flag == 0) ) // x(x-1)
      s = state[idx1] * (state[idx1] - 1) ;
    else // x*y or x*x
      s = state[idx1] * state[idx2] ;
  }

  return s ;
}

/* 
 *
 * compute the propensity ai of a reaction channel based on basis functions
 *
 * input: 
 *   channel_idx : 	index of the reaction channel
 *   state :		current state of the system 
 *   coeff_vec: 	current values of the unkonwn coefficients
 *
 * output: 
 *   s :		linear combination of the basis functions
 *
 */
double val_ai(int channel_idx, vector<int> &state, vector<vector<double> > & coeff_vec )
{
  int idx ;
  double s ; 
  s = 0 ; 
  // loop for each basis function used for this channel
  for ( int i = 0 ; i < coeff_vec[channel_idx].size() ; i ++ )
  {
    // index of basis function
    idx = basis_index_per_channel[channel_idx][i] ;
    // linear combination
    s += val_basis_funct(idx, state) * coeff_vec[channel_idx][i] ;
  }
  return s ;
}

/* 
 * print debug information
 *
 */
void dump_info(int idx, double tmp_ai, vector<int> & c_state, vector<vector<double> > & coeff_vec )
{
  int ii ;
  printf("\t State : ");
  for (int ii = 0 ; ii < n ; ii ++)
    printf("%d ", c_state[ii]); 
  printf("\n\t Channel idx: %d\n", idx) ;
  printf("\t Omega coefficients: ");
  for (int ii = 0 ; ii < coeff_vec[idx].size(); ii ++)
    printf("%.4lf\t", coeff_vec[idx][ii]); 
  printf("\n\t ai = %.4lf\n", tmp_ai);
}

/* 
 *
 * Compute (partial) log-likelihood function corresponding to certain reaction
 * channel
 *
 * If reaction types are known (know_reactions_flag=1), then function G(x)=x;
 * Otherwise, G(x) is a smooth function approximating \max(x,0).
 *
 * This function requires all processors to work together. 
 * Each processor does part of the computation using part of trajectory data.
 * The final result is obtained by calling MPI_Allreduce. 
 *
 * input:
 *
 *   i0 :		index of the reaction channel
 *   coeff_vec :	vector of coefficients  
 *
 * output:
 *   (partial) logarithmic of the likelihood function (divided by the total length
 *   of time)
 *
 */
double log_likelihood_partial( int i0, vector<vector<double> > & coeff_vec )
{
  double s, tmp_ai, ai, local_s ;
  int idx ;

  local_s = 0 ;
  // process trajectories that are distributed on the local processor
  for (int traj_idx = 0; traj_idx < local_N_traj; traj_idx ++)
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++) // loop for each jump (or reaction)
  {
    // subtract by 1, because index starts from one   
    idx = channel_idx_in_traj[traj_idx][i] - 1 ;

    // only when the jump is in channel i0 
    if (idx != i0) continue ;

    tmp_ai = val_ai(idx, traj_vec[traj_idx][i], coeff_vec) ;

    // since the channel occurs in trajectory, the propensity ai should be
    // positive
    if ( is_zero(tmp_ai) )
    {
      printf( "Error: ai=0! i=%d, idx=%d", i, idx ) ;

      if (mpi_rank == 0)
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
      tmp_ai = val_ai( i0, traj_vec[traj_idx][i], coeff_vec ) ;

    // compute the second part of the log-likelihood function
      local_s -= waiting_time_vec[traj_idx][i] * G(tmp_ai) ;
    }

  // sum up all processors
  MPI_Allreduce(&local_s, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;

  return s / total_T ;
}

/* 
 *
 * Compute gradient of the (partial) log-likelihood function corresponding a
 * certain reaction channel.
 *
 * Again, this function requires all processors to work together. 
 *
 * input:
 *   i0 	: 	index of the reaction channel
 *   coeff_vec  :	vector of coefficients  
 *
 * return:
 *   grad_coeff :	contains the gradient of the (partial) log-likelihood function (divided by the total length
 *   			of time)
 *
 */

void grad_log_likelihood_partial( int i0, vector<vector<double> > & coeff_vec, vector<vector<double> > & grad_coeff )
{
  double s, tmp_ai, s1, tmp , local_s ;
  int idx , basis_idx ;

  // compute each component of the gradient in the reaction channel i0
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

/* Shrinkage operator: 
*   	max(|x|-lambda, 0) * sgn(x)
*/
double shrinkage(double x, double lambda)
{
  if (x > lambda) 
    return x - lambda ;
  if (x < -lambda)
    return x + lambda ;
  return 0.0;
}

void print_grad_partial( int i, vector<vector<double> > & coeff_vec )
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

/* 
 * check whether a directory exists or not 
 *
 * input:
 *   dir_name :	path of the directory
 *
 * output:
 *   0, 	if the directory exists or is created
 *   -1,	if the directory doesn't exist and can not be created 
 *
 */
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

void p_l(int i, double L, vector<vector<double> > & yk, vector<vector<double> > & grad_f, vector<vector<double> > & vec_tmp) 
{
  double w , tmp, d1, d2, x, tol, residual ;

  for (int j = 0 ; j <  basis_index_per_channel[i].size() ; j ++)
  {
    w = regular_lambda * weights[i][j] ; 
    tmp = yk[i][j] - y grad_f[i][j] / L ;

    if ( is_zero(w) ) // if the weight is zero 
    {
      vec_tmp[i][j] = tmp ;
      continue ;
    }

    if (epsL1_flag == 0) // weighted l^1 norm
      vec_tmp[i][j] = shrinkage( tmp, w / L ) ;
    else  
    {  
      /* epsL1 : g(x) is w\times sqrt(x^2+\eps)
       *
       * Use newton method to solve: g(x) + 1/(2L) |x-(y-grad_f/L)|^2
       *
       * 1st derivative : wx/sqrt(x^2+\eps) + (x-(y-grad_f/L)/L
       * 2st derivative : 1/L + w\eps (x^2+\eps)^{-3/2}
       *
       * Starting from x0 = y-grad_f/L
       */

      x = tmp ;
      tol = 1e-10 ;
      residual = 1.0 ;
      while (residual > tol) 
      {
	d1 = w * x / sqrt(x*x+eps) + (x - tmp) / L ;
	d2 = 1.0 / L + w * eps / pow(x*x+eps, 1.5) ;
	// newton's method 
	x = x - d1 / d2 ;
	residual = fabs(d1 / d2) ;
      }

      if ((residual > tol) && (mpi_rank == 0))
      {
	 printf( "Warning: Newton's method doesn't converge!\n Residual=%.4e, Tolerance=%.4e\n", residual, tol) ;
	 fprintf( log_file, "Warning: Newton's method doesn't converge!\n Residual=%.4e, Tolerance=%.4e\n", residual, tol) ;
      }

      vec_tmp[i][j] = x ;
    }
  }
}


#include "sparse_learning.h"

/* 
 *
 * Determine channel information of each reaction from trajectory data.
 *
 *  This function is similar to the function:
 *
 *   	void find_channels_in_traj(vector<vector<vector<int> > > & traj_data)
 *
 *  in the file prepare.cpp. 
 *
 *  Different from that one in the file prepare.cpp, here: 
 *
 *  	1. both the list and the number of channel vectors are read from the file generated by ./prepare 
 *  	2. trajectories are distributed among processors
 *
 */

void determine_channel_index_of_each_reaction_in_traj(vector<vector<vector<int> > > & traj_data)
{
  map<vector<int>, int> channel_to_idx ;
  vector<int> vec_change ;
  int idx ;

  vec_change.resize(n) ;
  channel_idx_in_traj.resize(local_N_traj) ;
  Mi_in_traj.resize(local_N_traj) ;

  // indices of channels are indexed from 1 
  // channel_list is loaded from file, which is generated by running ./prepare
  for (int i = 0 ; i < channel_num ; i ++)
    channel_to_idx[channel_list[i]] = i + 1 ;

  // loop for each trajectory
  for (int traj_idx = 0 ; traj_idx < local_N_traj ; traj_idx ++)
  {
    channel_idx_in_traj[traj_idx].resize( num_state_in_traj[traj_idx] - 1 ) ;
    Mi_in_traj[traj_idx].resize(channel_num, 0) ;

    // for each reaction in the trajectory
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++)
    {
      // change vector
      for (int j = 0; j < n ; j ++)
	vec_change[j] = traj_data[traj_idx][i+1][j] - traj_data[traj_idx][i][j] ;

      if (channel_to_idx.count(vec_change) == 1)
      {
	idx = channel_to_idx[vec_change] ;
	channel_idx_in_traj[traj_idx][i] = idx ;
	// subtract by 1, because channel index starts from one   
	Mi_in_traj[traj_idx][idx - 1] ++ ;
      } else // something wrong...
      {
	printf("Error: reaction doesn't belong to any channels\n") ; 
	printf("\tchange vector: ");
	for (int ii = 0 ; ii < n ; ii++)
	{
	  printf( "%d ", vec_change[ii] );
	}
	printf("\n") ;
	exit(1) ;
      }
    }
  }

  int * l_mi_traj, * mi_traj_vec ;

  l_mi_traj = (int *) malloc( sizeof(int) * channel_num ) ;
  mi_traj_vec = (int *) malloc( sizeof(int) * channel_num ) ;

  for (int i = 0 ; i < channel_num; i ++)
    l_mi_traj[i] = 0 ;

  for (int i = 0 ; i < channel_num; i ++)
    for (int traj_idx = 0 ; traj_idx < local_N_traj ; traj_idx ++)
      l_mi_traj[i] += Mi_in_traj[traj_idx][i] ;
  MPI_Allreduce(l_mi_traj, mi_traj_vec, channel_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;

  Mi_in_all_traj.resize(channel_num) ;
  for (int i = 0 ; i < channel_num; i ++)
    Mi_in_all_traj[i] = mi_traj_vec[i] ;

  free(l_mi_traj) ;
  free(mi_traj_vec) ;

  // assign reactions to channels
  if (know_reactions_flag == 1)
  {
    reactions_in_channel.resize(channel_num) ;
    for (int i = 0; i < R; i ++)
    {
      if (channel_to_idx.count(vvec[i]) == 1)
      {
	idx = channel_to_idx[vvec[i]] ;
	reactions_in_channel[idx-1].push_back(i) ;
      }
      else 
      {
	if (mpi_rank == 0)
	{
	  printf("\nWarning: reaction %d doesn't correspond to any channel found in trajectories\n", i);
	  fprintf(log_file, "\nWarning: reaction %d doesn't correspond to any channel found in trajectories\n", i);
	}
      }
    }
  }
}

/* 
 * read the channel information from file
 *
 * Both the number of channels (channel_num) and the list of channel vectors (channel_vec) 
 * are read from the file 
 *
 * The user should run 
 *  	./prepare 
 * before run 
 * 	./sparse_learning 
 *
 */
void read_channels_info_from_file()
{
  char buf[100] ;
  int itmp ;
  vector<int> c_state ;
  ifstream in_file ;

  // this file is generated by running ./prepare first!
  sprintf( buf, "./output/channel_info.txt" ) ;
  if ( mpi_rank == 0 )
  {
    printf("\nReading channel information from file : %s\n\n", buf) ;
    fprintf(log_file, "\nReading channel information from file : %s\n\n", buf) ;
  }

  in_file.open(buf) ;

  // make sure the file is successfully open
  if ( (mpi_rank == 0) && (in_file.is_open() == 0) ) 
    {
      printf("Error: can not open file : %s.\n\tCommand ./prepare should be run first! \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s.\n\tCommand ./prepare should be run first! \n\n", buf) ;
      exit(1) ;
    }

  in_file >> channel_num >> itmp ;
  assert (itmp == n ) ;

  c_state.resize(n) ;
  channel_list.resize(0) ;

  // read the list of channels 
  for (int i = 0 ; i < channel_num ; i ++)
  {
    for (int j = 0 ; j < n ; j ++)
      in_file >>  c_state[j] ;
    channel_list.push_back(c_state) ;
  }

  if ( mpi_rank == 0 ) //  print the channel information
  {
    printf("No. of channels : %d\n", channel_num) ;
    fprintf(log_file, "No. of channels : %d\n", channel_num) ;

    for (int i =0 ; i < channel_num ; i ++)
    {
      printf("Change vector of the %dth channel :   ", i+1) ; 
      fprintf(log_file, "Change vector of the %dth channel :   ", i+1) ; 
      for (int j = 0 ; j < n ; j ++)
      {
	cout << channel_list[i][j] << ' ' ;
	fprintf( log_file, "%d ", channel_list[i][j] ) ;
      }
      cout << endl ;
      fprintf( log_file, "\n") ;
    }
  }

  in_file.close() ;
}

/*
 * Read trajectory data from file.
 * Each processor reads part of the whole data.
 *
 */
void read_trajectory_data() 
{
  char buf[100] ;
  ifstream in_file ;
  double t_now , tau , previous_tau , local_T ;
  int channel_idx ;
  string line ;
  vector<int> c_state ;

  num_state_in_traj.resize(local_N_traj) ;
  traj_vec.resize(local_N_traj) ;
  waiting_time_vec.resize(local_N_traj) ;
  t_vec.resize(local_N_traj) ;
  T_traj_vec.resize(local_N_traj) ;
  local_T = 0 ;

  // each process loads part of the trajectory data that it will process
  for ( int j = 0 ; j < local_N_traj ; j ++ )
    {
      sprintf( buf, "./traj_data/traj_%d.txt", j + local_traj_start_idx ) ;
      in_file.open(buf) ;
      if ( ! in_file.is_open() )
      {
	printf("Error: can not open file : %s, trajectory data is incomplete! \n\n", buf) ;
	exit(1) ;
      }

      in_file >> n ;
      c_state.resize(n) ;

      num_state_in_traj[j] = 0 ;
      // read line by line
      while ( getline(in_file, line) )
      {
	if (line.find_first_not_of(' ') != string::npos) 
	{
	  istringstream iss(line) ;
	  iss >> t_now ;
	  for (int i = 0 ; i < n ; i ++)
	    iss >> c_state[i] ;
	  iss >> tau ;

	  traj_vec[j].push_back(c_state) ;
	  t_vec[j].push_back(t_now) ;
	  waiting_time_vec[j].push_back(tau) ;

	  num_state_in_traj[j] ++ ;
	}
      }

      // the length of the current trajectory 
      T_traj_vec[j] = t_vec[j][num_state_in_traj[j] - 1] + waiting_time_vec[j][num_state_in_traj[j] - 1] ;
      // the total length of trajectories on the local processor
      local_T += T_traj_vec[j] ;

      in_file.close() ;
    }

  // compute the total length of all trajectories 
  MPI_Allreduce( &local_T, &total_T, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

  read_channels_info_from_file() ;

  determine_channel_index_of_each_reaction_in_traj(traj_vec) ;

  if (mpi_rank == 0)
  {
    printf( "Occurrence of each reaction channels within %d trajectories : [%d", N_traj, Mi_in_all_traj[0]) ;
    fprintf( log_file, "Occurrence of each reaction channels within %d trajectories : [%d", N_traj, Mi_in_all_traj[0]) ;
    for (int i = 1 ; i < channel_num ; i ++)
    {
      printf( ", %d", Mi_in_all_traj[i] ) ;
      fprintf( log_file, ", %d", Mi_in_all_traj[i] ) ;
    }

    printf( "]\n\nTotal length of time of %d trajectories : %.4f\n\n", N_traj, total_T ) ;
    fprintf( log_file, "]\n\nTotal length of time of %d trajectories : %.4f\n\n", N_traj, total_T ) ;
  }
}

/* 
 *
 * Read basis functions from file :    	
 * 	./output/basis_funct_info.txt
 *
 * The basis functions used for each channel will be read.
 *
 * The file is generated by running ./prepare
 *
 */
void read_basis_functions() 
{
  ifstream in_file ;
  int itmp ;
  char buf[100] ;

  // read basis functions from file, which is generated by running ./prepare!
  sprintf( buf, "./output/basis_funct_info.txt" ) ;
  in_file.open(buf) ;

  in_file >> itmp >> num_basis ;

  if ( (mpi_rank == 0) && (in_file.is_open() == 0) ) 
    {
      printf("Error: can not open file : %s.\n\tCommand ./prepare should be run first! \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s.\n\tCommand ./prepare should be run first! \n\n", buf) ;
      exit(1) ;
    }

  assert(itmp == n) ;

  // read basis functions
  basis_vec.resize(num_basis) ;
  for (int i = 0 ; i < num_basis ; i ++)
  {
    in_file >> itmp ;
    basis_vec[i].resize(itmp) ;
    for (int j = 0 ; j < itmp ; j ++)
      in_file >> basis_vec[i][j] ;
  }

  if (mpi_rank == 0)
  {
    printf("\n Loaded basis functions : %d\n", num_basis) ;
    //  read basis indices for each channel
    printf(" Reading basis functions for each channel...\n") ;

    fprintf(log_file, "\n Loaded basis functions : %d\n", num_basis) ;
    fprintf(log_file, " Reading basis functions for each channel...\n") ;
  }
  in_file >> itmp ;
  assert(itmp == channel_num) ;
  basis_index_per_channel.resize( channel_num ) ;

  total_unknown_omega_parameters = 0 ;
  // read number of basis used for each channel
  for ( int i = 0 ; i < channel_num ; i ++ )
  {
    in_file >> itmp ;
    assert(itmp <= num_basis) ;
    basis_index_per_channel[i].resize(itmp) ;
    total_unknown_omega_parameters += itmp ;
  }

  // read indices of basis functions for each channel
  for ( int i = 0 ; i < channel_num ; i ++ )
  {
    for ( int j = 0 ; j < basis_index_per_channel[i].size() ; j ++ )
    {
      in_file >> basis_index_per_channel[i][j] ;
    }
  }

  omega_vec.resize( channel_num ) ;
  omega_weights.resize( channel_num ) ;

  //  read sparsity weights of parameters for each channel
  if (mpi_rank == 0) 
  {
    printf(" Number of unknown parameters : %d\n", total_unknown_omega_parameters ) ;
    fprintf(log_file, " number of unknown parameters : %d\n", total_unknown_omega_parameters ) ;

    printf(" Reading (relative) weights of each parameter...\n") ;
    fprintf(log_file, " Reading (relative) weights of each parameter...\n") ;
  }

  for ( int i = 0 ; i < channel_num ; i ++ )
  {
    itmp = basis_index_per_channel[i].size() ; 
    omega_weights[i].resize(itmp) ;
    for (int j = 0 ; j < itmp ; j ++)
    {
      in_file >> omega_weights[i][j] ;
    }
  }

  //  read initial values of parameters for each channel
  if (mpi_rank == 0) 
  {
    printf(" Reading initial guesses ...\n") ;
    fprintf(log_file, " Reading initial guesses ...\n") ;
  }

  for ( int i = 0 ; i < channel_num ; i ++ )
  {
    itmp = basis_index_per_channel[i].size() ; 
    omega_vec[i].resize(itmp) ;
    for (int j = 0 ; j < itmp ; j ++)
      in_file >> omega_vec[i][j] ;
  }
}

/*
 * Output the solution of the optimization problem.
 *
 */

void output_omega()
{
  ofstream out_file ;
  char buf[100];
  sprintf( buf, "./output/omega_vec.txt" ) ;

  out_file.open(buf) ;
  if ( out_file.is_open() == 0 ) 
    {
      printf("Error: can not open file : %s. \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
      exit(1) ;
    }

  out_file << n << ' ' << channel_num << ' ' << num_basis << endl ;

  for (int i = 0; i < channel_num; i ++)
    out_file << omega_vec[i].size() << ' ';
  out_file << endl ;

  for (int i = 0; i < channel_num; i ++)
    for (int j = 0 ; j < omega_vec[i].size(); j ++)
      out_file << omega_vec[i][j] << ' ';

  out_file << endl ;
  out_file.close() ;
}

/*
 *
 * Scheme 1: iterative shrinkage-thresholding algorithm for solving the unknown parameters
 *
 * This iterative algorithm will be used if 
 * 	iter_scheme=0
 *
 * At each iteration step, the algorithm performs two updates:
 *
 *   1. update the values of parameters according to gradient descent
 *
 *   2. thresholding the updated parameters using the function
 * 		double shrinkage(double x, double lambda)
 *      in the file utils.cpp
 *
 */
void ista()
{
  double residual, tmp, tmp1 ;
  double new_omega ;
  int iter_step ; 
  vector<vector<double> > omega_grad_vec ;
  ofstream out_file ;
  char buf[100] ;

  // used in FISTA 
  vector<vector<double> > vec_tmp, previous_vec ;
  double L0, t1, eta, Lbar, t_new, t_old ;

  if (mpi_rank == 0)
  {
    sprintf( buf, "./output/iteration_omega_vec.txt" ) ;

    out_file.open(buf) ;
    if ( out_file.is_open() == 0 ) 
      {
	printf("Error: can not open file : %s. \n\n", buf) ;
	fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
	exit(1) ;
      }
    out_file << total_unknown_omega_parameters << endl ;
  }

  omega_grad_vec.resize( channel_num ) ;
  vec_tmp.resize( channel_num ) ;
  previous_vec.resize( channel_num ) ;
  for (int i = 0; i < channel_num; i ++)
  {
    vec_tmp[i].resize( basis_index_per_channel[i].size() ) ;
    yk[i].resize( basis_index_per_channel[i].size() ) ;
    omega_grad_vec[i].resize( basis_index_per_channel[i].size() ) ;
  }

  // initialize constants in FISTA 
  L0 = 1.0 ;
  t1 = 1.0 ;
  eta = 2.0 ;

  // solving the unknown coefficients for each channel
  for (int i =0 ; i < channel_num; i ++)
  {
    residual = 100.0 ;
    iter_step = 0 ;

    // initialize 
    t_old = t1 ;
    for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
      yk[i][j] = omega_vec[i][j] ;

    // update the parameters iteratively
    while ( (residual > stop_eps) && (iter_step < tot_step) ) 
    {
      residual = 0.0 ;
      Lbar = L0 ;

      // compute gradient of the log-likelihood functions
      grad_log_likelihood_partial(i, yk, omega_grad_vec) ;

      // evaluate the function at old point yk
      fval_old = -log_likelihood_partial(i, yk) ; 

      while (1) 
      {
	// projection by shrinkage
	p_l(i, Lbar, yk, omega_grad_vec, vec_tmp) ;

	// evaluate the function at new point
	fval_new = -log_likelihood_partial(i, vec_tmp) ; 

	// compute the function Q_L(x,y) ( without the term g(x)! ) 
	tmp = fval_old ;
	for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	  tmp += (vec_tmp[i][j] - yk[i][j]) * omega_grad_vec[i][j] + Lbar * 0.5 * (vec_tmp[i][j] - yk[i][j]) * (vec_tmp[i][j] - yk[i][j]) ;

	// check whether the condition is satisfied
	if (fval_new <= tmp) break ;

	// if not, increase the constant Lbar 
	Lbar *= eta ;
      }

      // update t_{k+1}
      t_new = (1 + sqrt(1 + 4 * t_old * t_old)) * 0.5 ;

      // update the vector yk 
      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	yk[i][j] = vec_tmp[i][j] + (t_old - 1) / t_new * ( vec_tmp[i][j] - omega_vec[i][j] ) ;

      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
      {
	// compute residual
	tmp = fabs( omega_vec[i][j] - vec_tmp[i][j] ) ;
	if (tmp > residual) residual = tmp ;
	  // update xk
	omega_vec[i][j] = vec_tmp[i][j] ;
      }

      iter_step ++ ;
      if (iter_step % output_interval == 0) // print information 
      {
	if (mpi_rank == 0)
	{
	  printf( "\nChannel idx=%d\tIteration step: %d\t \tResidual=%.6e ", i, iter_step, residual ) ;
	  fprintf( log_file, "\nChannel idx=%d\tIteration step: %d\t \tResidual=%.6e ", i, iter_step, residual ) ;
	}

	// all processors need to work together in order to compute log-likelihood
	tmp = -log_likelihood_partial(i, omega_vec) ; 
	tmp1 = l1_norm_partial(i,omega_vec, omega_weights) ;

	if (mpi_rank == 0)
	{
	  printf( "\nminus-likelihood = %.8f\t l1-norm : %.4f\t Cost = %.8f\n", tmp, tmp1, tmp + regular_lambda * tmp1 ) ;
	  fprintf( log_file, "\nminus-likelihood = %.8f\t l1-norm : %.4f\t Cost = %.8f\n", tmp, tmp1, tmp + regular_lambda * tmp1 ) ;
	  print_grad_partial(i, omega_vec) ;

	  out_file << "Channel idx=" << i << "\tIter_step=" << iter_step << endl ;
	  for (int j = 0 ; j < omega_vec[i].size(); j ++)
	    out_file << omega_vec[i][j] << ' ' ;
	  out_file << residual << endl ;
	}
      }
    }
  }

  if (mpi_rank == 0) 
  {
    output_omega() ;
    out_file.close() ;
  }
}

/* 
 * Scheme 2: iterative algorithm for solving the unknown parameters
 *
 * This iterative algorithm will be used if 
 * 	iter_scheme=1
 *
 * Different from Scheme 1, here the l^1 penalty is replaced by \sum_i w_i \times sqrt(|x_i|^2 + \eps)
 *
 * A simple gradient descent method is applied to minimize the smoothed optimization problem
 *
 */

void epsL1()
{
  double residual , tmp , new_omega, tmp1 ;
  int iter_step ; 
  vector<vector<double> > omega_grad_vec ;

  residual = 100.0 ;
  iter_step = 0 ;

  omega_grad_vec.resize( channel_num ) ;
  for (int i = 0; i < channel_num; i ++)
  {
    omega_grad_vec[i].resize( basis_index_per_channel[i].size() ) ;
  }

  while ( (residual > stop_eps) && (iter_step < tot_step) ) 
  {
    residual = 0.0 ;
    grad_log_likelihood(omega_vec, omega_grad_vec) ;

    for (int i =0 ; i < channel_num; i ++)
      for (int j = 0 ; j < omega_vec[i].size() ; j ++)
      {
	// gradient of the smoothed (l^1) penalty
	tmp = omega_weights[i][j] * omega_vec[i][j] / sqrt(omega_vec[i][j] * omega_vec[i][j] + eps) ;
	// gradient update 
	new_omega = omega_vec[i][j] - descent_dt * (omega_grad_vec[i][j] + regular_lambda * tmp) ;

	tmp = fabs( omega_vec[i][j] - new_omega ) ;
	if (tmp > residual) residual = tmp ;
	omega_vec[i][j] = new_omega ;
      }

    iter_step ++ ;
    if (iter_step % output_interval == 0) 
    {
      if (mpi_rank == 0)
      {
	printf( "\nIteration step: %d,\t \tResidual=%.6e ", iter_step, residual ) ;
	fprintf( log_file, "\nIteration step: %d,\t \tResidual=%.6e ", iter_step, residual ) ;
      }
      tmp = -log_likelihood(omega_vec) ;
      tmp1 = epsL1_norm(omega_vec, omega_weights) ;
      if (mpi_rank == 0)
      {
	printf( "\nminus-likelihood = %.8f\t epsL1-norm : %.4f\t Cost = %.8f\n", tmp, tmp1, tmp + regular_lambda * tmp1 ) ;
	fprintf( log_file, "\nminus-likelihood = %.8f\t epsL1-norm : %.4f\t Cost = %.8f\n", tmp, tmp1, tmp + regular_lambda * tmp1 ) ;
	print_grad(omega_vec) ;
      }
    }
  }
}

/*
 *
 * Directly estimation of the rate constants for reaction channels which contain
 * only one reaction.
 *
 * When the reaction types are known (know_reactions_flag=1), we are
 * minimizing the log-likelihood function itself. 
 *
 * In this case, when certain reaction channel only contains 1 reaction, then its unknown 
 * rate constant can be estimated directly by counting.  
 *
 * The results can used to compare with those obtained from iterative methods.
 *
 */

void direct_compute_channel_with_single_reaction() 
{
  int idx;
  double local_s, s ;

  if (mpi_rank == 0) cout << endl << "Rates of the following channels (contain one reaction) can be direct computed: " << endl << endl ;

  for ( int i = 0 ; i < channel_num ; i ++ )
    if ( basis_index_per_channel[i].size()==1 )
    {
      idx = basis_index_per_channel[i][0] ;
      local_s = 0 ;
      for (int traj_idx = 0 ; traj_idx < local_N_traj ; traj_idx ++)
	for (int j = 0; j < num_state_in_traj[traj_idx] ; j ++)
	  local_s += waiting_time_vec[traj_idx][j] * val_basis_funct(idx, traj_vec[traj_idx][j]) ;

      MPI_Allreduce( &local_s, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;

      if (mpi_rank == 0) 
      {
	printf("\tRate of channel %d: %.4f\n", i, Mi_in_all_traj[i] * 1.0 / s) ;
	fprintf(log_file, "\tRate of channel %d: %.4f\n", i, Mi_in_all_traj[i] * 1.0 / s) ;
      }
    }
}

int main ( int argc, char * argv[] ) 
{
  MPI_Init(&argc, &argv) ;
  char buf[50] ;

  clock_t start , end ;

  sprintf(buf, "./log/sparse_infer.log") ;
  if ( init(buf) < 0 ) return -1 ;

  start = clock() ;

  read_trajectory_data() ;

  read_basis_functions() ;

  // solve the parameters using different schemes
  switch (iter_scheme) {
  case 0 :  // iterative shrinkage-thresholding algorithm
    ista() ;
    break ;
  case 1 :  // epsL1
    epsL1() ;
    break ;
  }

  // when reaction types are known and a channel has only one reaction, then
  // the parameter can be computed directly
  if (know_reactions_flag == 1)
    direct_compute_channel_with_single_reaction() ;

  end = clock() ;

  if (mpi_rank == 0)
  {
    printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;
    fprintf(log_file, "\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;
    fclose(log_file) ;
  }

  MPI_Finalize() ; 

  return 0; 
}
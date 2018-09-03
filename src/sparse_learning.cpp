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

  // channel_list is loaded from file, which is generated by running ./prepare
  for (int i = 0 ; i < channel_num ; i ++)
    channel_to_idx[ channel_list[i] ] = i ;

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
	Mi_in_traj[traj_idx][idx] ++ ;
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

#if USE_MPI == 1
  MPI_Allreduce(l_mi_traj, mi_traj_vec, channel_num, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
#else
  for (int i = 0 ; i < channel_num; i ++)
    mi_traj_vec[i] = l_mi_traj[i] ;
#endif

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
	reactions_in_channel[idx].push_back(i) ;
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
#if USE_MPI == 1
  MPI_Allreduce( &local_T, &total_T, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
#else
  total_T = local_T ;
#endif

  if ( mpi_rank == 0 )
  {
    printf("\n========================================================\n") ;
    fprintf(log_file, "\n========================================================\n") ;
    printf("Reading channel information from file : %s\n", buf) ;
    fprintf(log_file, "Reading channel information from file : %s\n", buf) ;
  }

  read_channels_info_from_file() ;

  determine_channel_index_of_each_reaction_in_traj(traj_vec) ;

  if (mpi_rank == 0)
  {
    printf("\nIn total:  %d reaction channels\n\n", channel_num) ;
    fprintf(log_file, "\nIn total:  %d reaction channels\n\n", channel_num) ;

    for (int i =0 ; i < channel_num ; i ++)
    {
      printf("Channel %d: \tChange vector : (", i+1) ; 
      fprintf(log_file, "Channel %d: \tChange vector : (", i+1) ; 
      for (int j = 0 ; j < n ; j ++)
      {
	printf("%d", channel_list[i][j]) ;
	fprintf( log_file, "%d", channel_list[i][j] );
	if (j < n-1) {
	  printf(",");
	  fprintf(log_file, ",");
	} else 
	{
	  printf("),\t\t");
	  fprintf(log_file, "),\t\t");
	}
      }
      printf("Occurrence : %d\n", Mi_in_all_traj[i]) ;
      fprintf(log_file, "Occurrence : %d\n", Mi_in_all_traj[i]) ;
    }

    printf( "\nTotal length of time of %d trajectories : %.4f\n", N_traj, total_T ) ;
    fprintf( log_file, "\nTotal length of time of %d trajectories : %.4f\n", N_traj, total_T ) ;

    printf("========================================================\n\n") ;
    fprintf(log_file, "========================================================\n\n") ;
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
    printf("\n========================================================\n") ;
    fprintf(log_file, "\n========================================================\n") ;
    printf("Reading basis info from the file: ./output/basis_funct_info.txt\n" ) ;
    fprintf(log_file, "Reading basis info from the file: ./output/basis_funct_info.txt\n" ) ;
    printf("\nLoaded basis functions : %d\n", num_basis) ;
    //  read basis indices for each channel
    printf("Reading basis functions for each channel...\n") ;

    fprintf(log_file, "\nLoaded basis functions : %d\n", num_basis) ;
    fprintf(log_file, "Reading basis functions for each channel...\n") ;
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
    printf("Number of unknown parameters : %d\n", total_unknown_omega_parameters ) ;
    fprintf(log_file, "Number of unknown parameters : %d\n", total_unknown_omega_parameters ) ;

    printf("Reading (relative) weights of each parameter...\n") ;
    fprintf(log_file, "Reading (relative) weights of each parameter...\n") ;
  }

  for ( int i = 0 ; i < channel_num ; i ++ )
  {
    itmp = basis_index_per_channel[i].size() ; 
    omega_weights[i].resize(itmp) ;
    for (int j = 0 ; j < itmp ; j ++)
    {
      in_file >> omega_weights[i][j] ;

      if ( omega_weights[i][j] < 0 ) 
	{
	  if (mpi_rank == 0)
	  {
	    printf("Warning: %dth weight in the %dth channel is negative %.4f. Changed to 1.0\n", j, i, omega_weights[i][j] );
	    fprintf(log_file, "Warning: %dth weight in the %dth channel is negative %.4f. Changed to 1.0\n", j, i, omega_weights[i][j] );
	  }

	  omega_weights[i][j] = 1.0 ;
	}

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
      {
	in_file >> omega_vec[i][j] ;

	if ( (is_nonpositive(omega_vec[i][j]) == 1) && (know_reactions_flag == 1) )
	  {
	    if (mpi_rank == 0)
	    {
	      printf("Warning:  initial value of the %dth parameter in the %dth channel is nonpositive: %.4f. Changed to 1.0\n", j, i, omega_vec[i][j] ) ;
	      fprintf(log_file, "Warning:  initial value of the %dth parameter in the %dth channel is nonpositive: %.4f. Changed to 1.0\n", j, i, omega_vec[i][j] ) ;
	    }
	    omega_vec[i][j] = 1.0 ;
	  }
      }
  }

  if (mpi_rank == 0)
  {
    printf("========================================================\n\n") ;
    fprintf(log_file, "========================================================\n\n") ;
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
  {
    for (int j = 0 ; j < omega_vec[i].size(); j ++)
      out_file << std::setprecision(8) << omega_vec[i][j] << ' ';
    out_file << endl ;
  }

  out_file << endl ;
  out_file.close() ;

  printf("Results are written in the file: %s\n\n", buf) ;
  fprintf(log_file, "Results are written in the file: %s\n\n", buf) ;
}

/*
 * 
 * Iterative Shrinkage-Thresholding Algorithm (with backtracking) for solving the unknown parameters
 *
 * This function implements the ''ISTA with backtracking'' method in the paper:
 *
 *  	A. Beck and M. Teboulle, "A fast iterative shrinkage-thresholding algorithm for 
 *   linear inverse problems",  SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183–202, 2009. 
 *
 *
 * This function is similar to the function FISTA_backtracking(). It can be used to compare different numerical schemes. 
 *
 */
void ISTA_backtracking()
{
  double residual, tmp ;
  int iter_step ; 
  vector<vector<double> > omega_grad_vec ;
  ofstream out_file ;
  char buf[100] ;

  // used in ISTA 
  vector<vector<double> > vec_tmp ;
  double L0, eta, Lbar ;
  double fval_old, fval_new ;

  omega_grad_vec.resize( channel_num ) ;
  vec_tmp.resize( channel_num ) ;
  for (int i = 0; i < channel_num; i ++)
  {
    vec_tmp[i].resize( basis_index_per_channel[i].size() ) ;
    omega_grad_vec[i].resize( basis_index_per_channel[i].size() ) ;
  }

  // initialize constants in ISTA
  L0 = 1.0 ;
  eta = 2.0 ;

  /* 
   *
   * The log-likelihood function ln L(w) can be written as 
   * 	ln L(w) = \sum_{i=1}^K ln L_i(w_i), 
   * i.e., unknown parameters belonging to different reaction channels are decoupled.
   *
   * Therefore, we solve the unknown coefficients for one channel after another.
   *
   */
  for (int i =0 ; i < channel_num; i ++)
  {

    if (mpi_rank == 0)
    {
      printf("Solving coefficients for the reaction channel %d...\n", i) ;
      fprintf(log_file, "Solving coefficients for the reaction channel %d...\n", i) ;

      sprintf( buf, "./output/iteration_omega_vec_for_channel_%d.txt", i) ;

      out_file.open(buf) ;
      if ( out_file.is_open() == 0 ) 
	{
	  printf("Error: can not open file : %s. \n\n", buf) ;
	  fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
	  exit(1) ;
	}
      // number of unknown coefficients in the current channel
      out_file << basis_index_per_channel[i].size() << endl ;
    }

    // initialize 
    iter_step = 0 ;

    // update the parameters iteratively
    while ( iter_step < tot_step ) 
    {
      residual = 0.0 ;
      Lbar = L0 ;

      // compute the gradient of the log-likelihood functions
      grad_minus_log_likelihood_partial(i, omega_vec, omega_grad_vec) ;

      // evaluate the function at the old point x_{k-1}
      fval_old = minus_log_likelihood_partial(i, omega_vec) ; 

      /* 
       * compute Lbar
       *
       * After iteration, vec_tmp contains the updated state p_L(x_{k-1})
       */
      while (1) 
      {
	// projection by shrinkage
	p_L(i, Lbar, omega_vec, omega_grad_vec, vec_tmp) ;

	// evaluate the function at new point
	fval_new = minus_log_likelihood_partial(i, vec_tmp) ; 

	// compute the function Q_L(x,y) ( without the term g(x)! ) 
	tmp = fval_old ;
	for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	  tmp += (vec_tmp[i][j] - omega_vec[i][j]) * omega_grad_vec[i][j] + Lbar * 0.5 * (vec_tmp[i][j] - omega_vec[i][j]) * (vec_tmp[i][j] - omega_vec[i][j]) ;

	// check whether the condition is satisfied
	if (fval_new <= tmp) break ;

	// if not, increase the constant Lbar 
	Lbar *= eta ;
      }

      // compute residual
      residual = difference_of_two_vectors(omega_vec[i], vec_tmp[i]) / basis_index_per_channel[i].size() ;

      // check if the stop criteria is achieved
      if (residual < stop_eps) break ;

      // update x_k
      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	omega_vec[i][j] = vec_tmp[i][j] ;

      iter_step ++ ;
      if (iter_step % output_interval == 0) // print information 
      {
	if (mpi_rank == 0)
	{
	  printf( "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	  fprintf( log_file, "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	}

	tmp = penalty_g_partial(i,omega_vec, omega_weights) ;

	if (mpi_rank == 0)
	{
	  printf( "\n\tminus-likelihood = %.8e\t\t penalty g_i(x) = %.5e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;
	  fprintf( log_file, "\n\tminus-likelihood = %.8e\t\t penalty g_i(x)= %.4e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;

	  print_omega_coefficients(i, omega_vec) ;

	  out_file << iter_step << "\t" << std::setprecision(8) ;

	  for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	    out_file << omega_vec[i][j] << ' ' ;

	  out_file << "\t" << std::setprecision(8) << fval_new + tmp << "\t" << residual << endl ;
	}
      }
    }

    if (mpi_rank == 0)
    {
      printf("\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;
      fprintf(log_file, "\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;

      out_file.close() ;
    }
  }

  tmp = 0 ;
  for (int i = 0 ; i < channel_num ; i ++)
  {
    tmp += minus_log_likelihood_partial( i, omega_vec ) ;
    tmp += penalty_g_partial( i, omega_vec, omega_weights ) ;
  }

  if (mpi_rank == 0) 
  {
    printf( "Final cost = %.6f\n\n", tmp ) ;
    fprintf( log_file, "Final cost = %.6f\n\n", tmp ) ;

    output_omega() ;
  }
}

/*
 * Fast Iterative Shrinkage-Thresholding Algorithm (with backtracking) for solving the unknown parameters
 *
 * This function implements the FISTA method introduced in the paper:
 *
 *  	A. Beck and M. Teboulle, "A fast iterative shrinkage-thresholding algorithm for 
 *   linear inverse problems",  SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183–202, 2009. 
 *
 */
void FISTA_backtracking()
{
  double residual, tmp ;
  int iter_step ; 
  vector<vector<double> > omega_grad_vec ;
  ofstream out_file ;
  char buf[100] ;

  // used in FISTA 
  vector<vector<double> > vec_tmp, yk ;
  double L0, t1, eta, Lbar, t_new, t_old ;
  double fval_old, fval_new ;

  omega_grad_vec.resize( channel_num ) ;
  vec_tmp.resize( channel_num ) ;
  yk.resize( channel_num ) ;
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

  /* 
   *
   * The log-likelihood function ln L(w) can be written as 
   * 	ln L(w) = \sum_{i=1}^K ln L_i(w_i), 
   * i.e., unknown parameters belonging to different reaction channels are decoupled.
   *
   * Therefore, we solve the unknown coefficients for one channel after another.
   *
   */
  for (int i =0 ; i < channel_num; i ++)
  {

    if (mpi_rank == 0)
    {
      printf("Solving coefficients for the reaction channel %d...\n", i) ;
      fprintf(log_file, "Solving coefficients for the reaction channel %d...\n", i) ;

      sprintf( buf, "./output/iteration_omega_vec_for_channel_%d.txt", i) ;

      out_file.open(buf) ;
      if ( out_file.is_open() == 0 ) 
	{
	  printf("Error: can not open file : %s. \n\n", buf) ;
	  fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
	  exit(1) ;
	}
      // number of unknown coefficients in the current channel
      out_file << basis_index_per_channel[i].size() << endl ;
    }

    // initialize 
    t_old = t1 ;
    iter_step = 0 ;

    for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
      yk[i][j] = omega_vec[i][j] ;

    // update the parameters iteratively
    while ( iter_step < tot_step ) 
    {
      residual = 0.0 ;
      Lbar = L0 ;

      // compute the gradient of the log-likelihood functions
      grad_minus_log_likelihood_partial(i, yk, omega_grad_vec) ;

      // evaluate the function at old point yk
      fval_old = minus_log_likelihood_partial(i, yk) ; 

      /* 
       * compute Lbar
       *
       * After iteration, vec_tmp contains the updated state p_L(yk)
       */
      while (1) 
      {
	// projection by shrinkage
	p_L(i, Lbar, yk, omega_grad_vec, vec_tmp) ;

	// evaluate the function at new point
	fval_new = minus_log_likelihood_partial(i, vec_tmp) ; 

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

      // compute residual
      residual = difference_of_two_vectors(omega_vec[i], vec_tmp[i]) / basis_index_per_channel[i].size() ;

      // check if the stop criteria is achieved
      if (residual < stop_eps) break ;

      // update x_k
      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	omega_vec[i][j] = vec_tmp[i][j] ;

      // update t_k
      t_old = t_new ;

      iter_step ++ ;
      if (iter_step % output_interval == 0) // print information 
      {
	if (mpi_rank == 0)
	{
	  printf( "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	  fprintf( log_file, "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	}

	tmp = penalty_g_partial(i,omega_vec, omega_weights) ;

	if (mpi_rank == 0)
	{
	  printf( "\n\tminus-likelihood = %.8e\t\t penalty g_i(x) = %.5e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;
	  fprintf( log_file, "\n\tminus-likelihood = %.8e\t\t penalty g_i(x)= %.4e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;

	  print_omega_coefficients(i, omega_vec) ;

	  out_file << iter_step << "\t" << std::setprecision(8) ;

	  for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	    out_file << omega_vec[i][j] << ' ' ;

	  out_file << "\t" << std::setprecision(8) << fval_new + tmp << "\t" << residual << endl ;
	}
      }
    }

    if (mpi_rank == 0)
    {
      printf("\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;
      fprintf(log_file, "\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;

      out_file.close() ;
    }
  }

  tmp = 0 ;
  for (int i = 0 ; i < channel_num ; i ++)
  {
    tmp += minus_log_likelihood_partial( i, omega_vec ) ;
    tmp += penalty_g_partial( i, omega_vec, omega_weights ) ;
  }

  if (mpi_rank == 0) 
  {
    printf( "Final cost = %.6f\n\n", tmp ) ;
    fprintf( log_file, "Final cost = %.6f\n\n", tmp ) ;

    output_omega() ;
  }
}

/*
 *
 * Simple gradient descent method (with fixed step size) to solve the smooth optimization problem,
 * i.e., when either epsL1_flag=1 or know_reactions_flag=1.
 * (in the latter case there is no penalty term)
 *
 */

void grad_descent_smooth() 
{
  double residual, tmp ;
  int iter_step ; 
  vector<vector<double> > omega_grad_vec ;
  ofstream out_file ;
  char buf[100] ;

  vector<vector<double> > vec_tmp ;
  double fval_old, fval_new ;

  omega_grad_vec.resize( channel_num ) ;
  vec_tmp.resize( channel_num ) ;
  for (int i = 0; i < channel_num; i ++)
  {
    vec_tmp[i].resize( basis_index_per_channel[i].size() ) ;
    omega_grad_vec[i].resize( basis_index_per_channel[i].size() ) ;
  }

  /* 
   *
   * The log-likelihood function ln L(w) can be written as 
   * 	ln L(w) = \sum_{i=1}^K ln L_i(w_i), 
   * i.e., unknown parameters belonging to different reaction channels are decoupled.
   *
   * Therefore, we solve the unknown coefficients for one channel after another.
   *
   */
  for (int i =0 ; i < channel_num; i ++)
  {

    if (mpi_rank == 0)
    {
      printf("Solving coefficients for the reaction channel %d...\n", i) ;
      fprintf(log_file, "Solving coefficients for the reaction channel %d...\n", i) ;

      sprintf( buf, "./output/iteration_omega_vec_for_channel_%d.txt", i) ;

      out_file.open(buf) ;
      if ( out_file.is_open() == 0 ) 
	{
	  printf("Error: can not open file : %s. \n\n", buf) ;
	  fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
	  exit(1) ;
	}
      // number of unknown coefficients in the current channel
      out_file << basis_index_per_channel[i].size() << endl ;
    }

    // initialize 
    iter_step = 0 ;

    // update the parameters iteratively
    while ( iter_step < tot_step ) 
    {
      residual = 0.0 ;

      // compute the gradient of the log-likelihood functions
      grad_minus_log_likelihood_partial(i, omega_vec, omega_grad_vec) ;

      // update the vector by gradient descent
      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
      {
	tmp = omega_vec[i][j] ;
	vec_tmp[i][j] = omega_vec[i][j] - grad_dt * ( omega_grad_vec[i][j] + regular_lambda * omega_weights[i][j] * tmp / sqrt(tmp * tmp + eps) ) ;
      }

      // compute residual
      residual = difference_of_two_vectors(omega_vec[i], vec_tmp[i]) / basis_index_per_channel[i].size() ;

      // check if the stop criteria is achieved
      if (residual < stop_eps) break ;

      // update x_k
      for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	omega_vec[i][j] = vec_tmp[i][j] ;

      iter_step ++ ;
      if (iter_step % output_interval == 0) // print information 
      {
	if (mpi_rank == 0)
	{
	  printf( "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	  fprintf( log_file, "\nChannel idx = %d\t\tIteration step = %d\t \tResidual = %.6e ", i, iter_step, residual ) ;
	}

	tmp = penalty_g_partial(i,omega_vec, omega_weights) ;

	if (mpi_rank == 0)
	{
	  printf( "\n\tminus-likelihood = %.8e\t\t penalty g_i(x) = %.5e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;
	  fprintf( log_file, "\n\tminus-likelihood = %.8e\t\t penalty g_i(x)= %.4e\t\t Cost = %.8e\n", fval_new, tmp, fval_new + tmp ) ;

	  print_omega_coefficients(i, omega_vec) ;

	  out_file << iter_step << "\t" << std::setprecision(8) ;

	  for (int j = 0 ; j < basis_index_per_channel[i].size() ; j ++)
	    out_file << omega_vec[i][j] << ' ' ;

	  out_file << "\t" << std::setprecision(8) << fval_new + tmp << "\t" << residual << endl ;
	}
      }
    }

    if (mpi_rank == 0)
    {
      printf("\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;
      fprintf(log_file, "\nSolving coefficients for the reaction channel %d... finished.\nTotal iteration steps=%d,\t Final residual=%.4e\n\n", i, iter_step, residual ) ;

      out_file.close() ;
    }
  }

  tmp = 0 ;
  for (int i = 0 ; i < channel_num ; i ++)
  {
    tmp += minus_log_likelihood_partial( i, omega_vec ) ;
    tmp += penalty_g_partial( i, omega_vec, omega_weights ) ;
  }

  if (mpi_rank == 0) 
  {
    printf( "Final cost = %.6f\n\n", tmp ) ;
    fprintf( log_file, "Final cost = %.6f\n\n", tmp ) ;

    output_omega() ;
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

  if (mpi_rank == 0) 
    {
      printf("\n========================================================\n") ;
      fprintf(log_file, "\n========================================================\n") ;

      printf("Direct calculation of rates of the following channels (which contains one reaction):\n") ;
      fprintf(log_file, "Direct calculation of rates of the following channels (which contains one reaction):\n") ;
    }

  for ( int i = 0 ; i < channel_num ; i ++ )
    if ( basis_index_per_channel[i].size()==1 )
    {
      idx = basis_index_per_channel[i][0] ;
      local_s = 0 ;
      for (int traj_idx = 0 ; traj_idx < local_N_traj ; traj_idx ++)
	for (int j = 0; j < num_state_in_traj[traj_idx] ; j ++)
	  local_s += waiting_time_vec[traj_idx][j] * val_basis_funct(idx, traj_vec[traj_idx][j]) ;

#if USE_MPI == 1
      MPI_Allreduce( &local_s, &s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) ;
#else
      s = local_s ;
#endif

      if (mpi_rank == 0) 
      {
	printf("\tRate of channel %d: %.8f\n", i, Mi_in_all_traj[i] * 1.0 / s) ;
	fprintf(log_file, "\tRate of channel %d: %.8f\n", i, Mi_in_all_traj[i] * 1.0 / s) ;
      }
    }
  if (mpi_rank == 0)
  {
    printf("========================================================\n") ;
    fprintf(log_file, "========================================================\n") ;
  }
}

void check_solver_id()
{
  if ( (solver_id > 3) || (solver_id < 1) )
    {
      if (mpi_rank == 0)
      {
	printf( "Warning: No such solver (solver_id=%d)!  Reset to solver_id=1 (FISTA method) \n\n", solver_id) ;
	fprintf( log_file, "Warning: No such solver (solver_id=%d)!  Reset to solver_id=1 (FISTA method) \n\n", solver_id) ;
      }
      solver_id = 1 ;
    }

  if (solver_id == 3) 
  {
    if ( (know_reactions_flag==0) && (epsL1_flag==0) )
    {
      if (mpi_rank == 0)
      {
	printf( "Warning:  The gradient descent method (solver_id=3) can only be used when the object function is smooth, not for l^1 sparsity optimization. Reset solver_id=1 (FISTA method) \n\n") ;
	fprintf( log_file, "Warning:  The gradient descent method (solver_id=3) can only be used when the object function is smooth, not for l^1 sparsity optimization. Reset solver_id=1 (FISTA method) \n\n") ;
      }
      solver_id = 1 ;
    }  
    else 
      if ( (is_nonpositive(grad_dt)==1) && (mpi_rank == 0) )
	  {
	    printf("Error: step-size in the gradient descent method=%.4f, It should be positve!\n\n", grad_dt) ; 
	    fprintf(log_file, "Error: step-size in the gradient descent method=%.4f, It should be positve!\n\n", grad_dt) ; 
	    exit(1) ;
	  }
  }

  // print information about the solver 
  if (mpi_rank == 0) 
  {
    printf("solver_id=%d\n", solver_id);
    fprintf(log_file, "solver_id=%d\n", solver_id);

    switch (solver_id) {
    case 1:
	printf( "Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) will be used to solve the coefficients.\n\n" ) ; 
	fprintf( log_file, "Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) will be used to solve the coefficients.\n\n" ) ; 
	break ;
    case 2:
	printf( "Iterative Shrinkage-Thresholding Algorithm (ISTA) will be used to solve the coefficients.\n\n" ) ; 
	fprintf( log_file, "Iterative Shrinkage-Thresholding Algorithm (ISTA) will be used to solve the coefficients.\n\n" ) ; 
	break ;
    case 3:
	printf("gradient descent method will be used. step-size=%.4f.\n\n", grad_dt) ; 
	fprintf( log_file, "gradient descent method will be used. step-size=%.4f.\n\n", grad_dt) ; 
	break ;
    }
  }
}

int main ( int argc, char * argv[] ) 
{

#if USE_MPI == 1
  MPI_Init(&argc, &argv) ;
#endif

  char buf[50] ;

  clock_t start , end ;

  start = clock() ;

  sprintf(buf, "./log/sparse_infer.log") ;
  if ( init(buf) < 0 ) return -1 ;

  check_solver_id() ;

  read_trajectory_data() ;

  read_basis_functions() ;

  switch (solver_id) {
    case 1 :
      // solve the parameters using the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA)
    	FISTA_backtracking() ;
	break; 
    case 2 :
      // solve the parameters using the Iterative Shrinkage-Thresholding Algorithm (ISTA)
  	ISTA_backtracking() ;
	break; 
    case 3 :
	/* 
	 * Simple (and slow) gradient descent method. 
	 *
	 * Only when either epsL1_flag=1 or know_reactions_flag=1, i.e., the object function is smooth.
	 *
	 */
	grad_descent_smooth() ;
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

#if USE_MPI == 1
  MPI_Finalize() ; 
#endif

  return 0; 
}

#include "sparse_learning.h"

/*
 * Determine all reaction channels from the N_traj trajectories data.
 *
 * The following data will be updated:
 *
 * channel_num:  		number of reaction channels
 *
 * channel_list: 		list of reaction channels (each channel is represented as a
 * 		 		vector (dimension n).
 *
 * channel_idx_in_traj: 	index of reaction channel for each jump in each
 * 			 	trajectory 
 *
 * Mi_in_traj: 			for every trajectory, number of occurrence of each reaction channel 
 * 	       			is recorded.
 *
 * Mi_in_all_traj: 		occurrence number of each reaction channel
 * 				within all trajectories 
 *
 */

void find_channels_in_traj(vector<vector<vector<int> > > & traj_data)
{
  map<vector<int>, int> channel_to_idx ;
  set<vector<int> > channel_vec_set ;
  vector<int> vec_change ;
  int idx ;

  vec_change.resize(n) ;
  channel_num = 0 ;
  channel_idx_in_traj.resize(N_traj) ;
  Mi_in_traj.resize(N_traj) ;

  // loop for each trajectory 
  for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
  {
    channel_idx_in_traj[traj_idx].resize( num_state_in_traj[traj_idx] - 1 ) ;

    // for each reaction in the trajectory 
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++)
    {
      // compute change vector 
      for (int j = 0; j < n; j ++)
	vec_change[j] = traj_data[traj_idx][i+1][j] - traj_data[traj_idx][i][j] ;

      // add a new channel to the set, if the change vector doesn't belong to the existing ones
      if ( channel_vec_set.count(vec_change) == 0 ) channel_vec_set.insert(vec_change) ;
    }
  }

  // channel index starts from 0
  // create list of channels 
  channel_num = 0 ;
  for (set<vector<int> >::iterator it = channel_vec_set.begin() ; it != channel_vec_set.end() ; it ++)
    {
      channel_to_idx[*it] = channel_num ;
      channel_list.push_back(*it) ;
      channel_num ++ ;
    }

  for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
  {
    Mi_in_traj[traj_idx].resize(channel_num, 0) ;
    // for each reaction in the trajectory 
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++)
    {
      // compute change vector
      for (int j = 0; j < n; j ++)
	vec_change[j] = traj_data[traj_idx][i+1][j] - traj_data[traj_idx][i][j] ;

      // determine the index for each reaction
      idx = channel_to_idx[vec_change] ;
      channel_idx_in_traj[traj_idx][i] = idx ;
      // increase counter 
      Mi_in_traj[traj_idx][idx] ++ ;
    }
  }

  Mi_in_all_traj.resize(channel_num, 0) ;
  // sum up all trajectories
  for (int i = 0 ; i < channel_num; i ++)
    for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
      Mi_in_all_traj[i] += Mi_in_traj[traj_idx][i] ;
}

/* 
 * Output the information of channels to file. This information will be
 * further used by ./sparse_learning. 
 * 
 */
void output_channel_info()
{
  ofstream out_file ;
  char buf[100] ;
  sprintf( buf, "./output/channel_info.txt") ;

  printf("\nInformation of channels are written into file: %s\n", buf) ;
  fprintf(log_file, "\nInformation of channels are written into file: %s\n", buf) ;

  out_file.open(buf) ;

  // make sure the file is open successfully
  if ( ! out_file.is_open() )
    {
      printf("Error: can not open file : %s \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s \n\n", buf) ;
      exit(1) ;
    }

  out_file << channel_num << ' ' << n << endl ;

  // output change vector of each channel
  for (int i = 0; i < channel_num; i ++)
  {
    for (int j = 0; j < n ; j++)
      out_file << channel_list[i][j] << ' ';
    out_file << endl ;
  }

  out_file << endl ;

  // for each trajectory, output the occurrence number of each channel
  for (int i = 0; i < N_traj; i ++)
  {
    for (int j = 0; j < channel_num ; j++)
    {
      out_file << Mi_in_traj[i][j] << ' ';
    }
    out_file << endl ;
  }

  out_file.close() ;
}

/*
 *
 * Process trajectory data under directory ./traj_data, which were generated by
 * ./ssa using SSA method.
 *
 * traj_vec : 		vector containing all trajectories  
 * t_vec : 		vector of time for each state in each trajectory
 * waiting_time_vec : 	vector of waiting time for each reaction in each trajectory 
 * num_state_in_traj :	number of total states in each trajectory 
 *
 */
void process_data() 
{
  char buf[100] ;
  ifstream in_file ;
  double t_now , tau , previous_tau ;
  int channel_idx ;
  string line ;
  vector<int> c_state ;

  num_state_in_traj.resize(N_traj) ;
  traj_vec.resize(N_traj) ;
  waiting_time_vec.resize(N_traj) ;
  t_vec.resize(N_traj) ;

  printf("Number of trajectories = %d\n\n", N_traj) ;
  fprintf(log_file , "Number of trajectories = %d\n\n", N_traj) ;

  tot_num_reaction_in_traj = 0 ;

  // read trajectory data from file 
  for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
  {
    sprintf( buf, "./traj_data/traj_%d.txt", traj_idx ) ;

    printf("processing %dth trajectory...\n", traj_idx) ;
    fprintf(log_file, "processing %dth trajectory...\n", traj_idx) ;

    in_file.open(buf) ;
    if ( ! in_file.is_open() )
    {
      printf("Error: can not open file : %s, trajectory data is incomplete! \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s, trajectory data is incomplete! \n\n", buf) ;
      exit(1) ;
    }

    in_file >> n ;
    c_state.resize(n) ;

    num_state_in_traj[traj_idx] = 0 ;

    // read the file line by line
    while ( getline(in_file, line) )
    {
      if (line.find_first_not_of(' ') != string::npos) 
      {
	istringstream iss(line) ;
	iss >> t_now ;
	for (int i = 0 ; i < n ; i ++)
	  iss >> c_state[i] ;
	iss >> tau ;

	traj_vec[traj_idx].push_back(c_state) ;
	t_vec[traj_idx].push_back(t_now) ;
	waiting_time_vec[traj_idx].push_back(tau) ;

	num_state_in_traj[traj_idx] ++ ;
      }
    }

    tot_num_reaction_in_traj += (num_state_in_traj[traj_idx] - 1) ;
    in_file.close() ;
  }

  // determine channels from trajectory data
  find_channels_in_traj(traj_vec) ;

  printf("\n\n========================================================\n") ;
  fprintf(log_file, "\n\n========================================================\n") ;

  printf("In total:\n   %d reactions,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;
  fprintf(log_file, "In total:\n   %d reactions,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;

  for (int i =0 ; i < channel_num ; i ++)
  {
    printf("Channel %d: \tChange vector : (", i) ; 
    fprintf(log_file, "Channel %d: \tChange vector : (", i) ; 
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

  output_channel_info() ;

  printf("========================================================\n\n") ;
  fprintf(log_file, "========================================================\n\n") ;
}

void find_min_max_val_of_basis_in_traj(int c_idx, vector<double> & min_val, vector<double> & max_val) 
{
  int idx ;
  double tmp ;

  for (int j = 0 ; j < num_basis ; j ++)
  {
    min_val[j] = 1e12 ;
    max_val[j] = -1e12 ;
    // loop for each trajectory 
    for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
      for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++) // loop for each jump (or reaction)
    {
      // get the index of the channel
      idx = channel_idx_in_traj[traj_idx][i] ;

      // only when the jump is in channel i0 
      if (idx != c_idx) continue ;

      // compute the value of basis function at the current state
      tmp = val_basis_funct(j, traj_vec[traj_idx][i]) ;

      if (tmp < min_val[j]) min_val[j] = tmp ;
      if (tmp > max_val[j]) max_val[j] = tmp ;
    }
  }
}

/* 
 *
 * Write basis functions to file. The information of basis functions will be
 * used by ./sparse_learning
 *
 * In the current implementation, only 1st/2nd order reactions are
 * considered. Correspondingly, the basis functions are either 1st/2nd
 * polynomials. 
 *
 * The number of basis functions (num_basis) depends on the value of poly_order. 
 * The task of the code ./sparse_learning is to approximate propensity functions of all
 * channels using basis functions. 
 *
 * There are n 1st order polynomial basis: 
 * 	x_1, x_2, ..., x_{n}
 *
 * There are n*(n+1)/2 2nd order polynomial basis, since 
 * 	x_1*x_1, x_1*x_2, ... x_{n-1}*x_n, x_n*x_n
 *
 * If poly_order=1, then only 1st order polynomial basis are used, i.e., num_basis=n;
 *
 * If poly_order=2, then both 1st and 2nd order polynomial basis are used, i.e., num_basis=n*(n+3)/2 ;
 *
 * By default, all basis functions will be used to approximate propensity
 * function of each reaction channel. 
 *
 * The (sparse) weights of rate constants are set to 1.0 in the output file.
 *
 */
void write_basis_functions() 
{
  vector<int> reactant_idx ;

  // initialize the vector before including any basis functions
  basis_vec.resize(0) ;

  // include constant function one as basis function
  reactant_idx.resize(0) ; 
  basis_vec.push_back( reactant_idx ) ;

  // include polynomial functions of order 1 (linear) as basis functions
  reactant_idx.resize(1) ; 
  for (int i =0 ; i < n; i++)
  {
    reactant_idx[0] = i ;
    basis_vec.push_back( reactant_idx ) ;
  }

  // include quadratic functions as basis functions
  if (poly_order == 2)
  {
    reactant_idx.resize(2) ; 
    for (int i=0 ; i < n ; i++)
      for (int j=i ; j < n ; j++)
      {
	reactant_idx[0] = i ;
	reactant_idx[1] = j ;
	basis_vec.push_back( reactant_idx ) ;
      }
    num_basis = 1 + n * (3 + n) / 2 ;   
  } else 
  {
    num_basis = 1 + n ; 
  }

  // initialization
  omega_vec.resize(channel_num) ;
  for (int i = 0 ; i < channel_num ; i ++)
    omega_vec[i].resize( num_basis, 0.0 ) ;

  printf("\n\n========================================================\n") ;
  fprintf(log_file, "\n\n========================================================\n") ;

  // print basis functions into file
  printf("Number of basis functions = %d\n", num_basis) ;
  fprintf(log_file, "Number of basis functions = %d\n", num_basis) ;

  ofstream out_file ;
  char buf[100] ;
  sprintf( buf, "./output/basis_funct_info.txt" ) ;

  printf("\nBasis functions are written into the file: %s\n\n", buf) ;
  fprintf(log_file, "\nBasis functions are written into the file: %s\n\n", buf) ;
  out_file.open(buf) ;
  if ( ! out_file.is_open() )
    {
      printf("Error: can not open output file : %s. \n\n", buf) ;
      fprintf(log_file, "Error: can not open output file : %s. \n\n", buf) ;
      exit(1) ;
    }

  out_file << n << ' ' << num_basis << endl ;
  for (int i = 0 ; i < num_basis ; i ++)
  {
    out_file << basis_vec[i].size() << ' ' ;
    for (int j = 0 ; j < basis_vec[i].size() ; j ++)
      out_file << basis_vec[i][j] << ' ' ;
    out_file << endl ;
  }
  out_file << endl ;
  out_file.close() ;

  vector<double> min_val , max_val ;

  min_val.resize( num_basis ) ;
  max_val.resize( num_basis ) ;

  // write basis information of each channel 
  for (int i = 0 ; i < channel_num ; i ++)
  {
    sprintf( buf, "./output/basis_of_channel_%d.txt", i) ;
    printf("Basis info of channel %d is written to : %s\n", i, buf) ;
    fprintf(log_file, "Basis info of channel %d is written to : %s\n", i, buf) ;
    out_file.open(buf) ;
    if ( ! out_file.is_open() )
      {
	printf("Error: can not open output file : %s. \n\n", buf) ;
	fprintf(log_file, "Error: can not open output file : %s. \n\n", buf) ;
	exit(1) ;
      }

    find_min_max_val_of_basis_in_traj(i, min_val, max_val) ;

    // output number of basis functions for channel i, 
    // all basis functions will be used by default
    out_file << num_basis << ' ' << endl ;

    for (int j = 0 ; j < num_basis ; j ++)
      out_file << j << ' ' ;
    out_file << endl ;

    // output weight of each parameter 
    for (int j = 0 ; j < num_basis ; j ++)
      out_file << 1.0 << ' ' ;
    out_file << endl ;

    // initial value of parameters for each basis functions 
    for (int j = 0 ; j < num_basis ; j ++)
      out_file << 0.0 << ' ' ;
    out_file << endl << endl ;

    // min value of each basis function 
    for (int j = 0 ; j < num_basis ; j ++)
      out_file << min_val[j] << ' ' ;
    out_file << endl ;

    // max value of each basis function 
    for (int j = 0 ; j < num_basis ; j ++)
      out_file << max_val[j] << ' ' ;
    out_file << endl ;

    out_file.close() ;

    printf("Range of each basis function in channel %d: \n", i) ;
    fprintf(log_file, "Range of each basis function in channel %d: \n", i) ;
    for (int j = 0 ; j < num_basis ; j ++)
    {
      printf("[%.2e, %.2e]\t", min_val[j], max_val[j]) ;
      fprintf(log_file, "[%.2e, %.2e]\t", min_val[j], max_val[j]) ;
    }

    printf("\n\n");
    fprintf(log_file, "\n\n");
  }

  /* 
   * Only indices contained in the following file will be learned.
   * By default, indices of all channels are included.
   */
  sprintf( buf, "./output/channels_to_learn.txt" ) ;
  out_file.open(buf) ;
  if ( ! out_file.is_open() )
    {
      printf("Error: can not open output file : %s. \n\n", buf) ;
      fprintf(log_file, "Error: can not open output file : %s. \n\n", buf) ;
      exit(1) ;
    }

  out_file << channel_num << endl ;

  for (int i = 0 ; i < channel_num; i ++)
    out_file << i << ' ';

  out_file << endl ;
  out_file.close() ;

  printf("\nIndices of channels which will be learned are written to: %s\n", buf);
  fprintf(log_file, "\nIndices of channels which will be learned are written to: %s\n", buf);

  printf("========================================================\n\n") ;
  fprintf(log_file, "========================================================\n\n") ;

}

int main ( int argc, char * argv[] ) 
{
#if USE_MPI == 1
  MPI_Init(&argc, &argv) ;
#endif

  clock_t start , end ;
  char buf[30]; 

  know_reactions_flag = 0 ;

  sprintf(buf, "./log/prepare.log") ;
  if (init(buf) < 0) return -1 ;

  // this program can only run sequentially
  assert(mpi_size == 1) ;

  // make sure the directory is ready
  sprintf(buf, "./output") ;
  if (dir_check(buf) != 0) 
    {
      printf ("\nError: directory %s can not be created. Fix the error manually!\n\n", buf) ;
      fprintf (log_file, "\nError: directory %s can not be created. Fix the error manually!\n\n", buf) ;
      exit(1) ;
    }

  start = clock() ;

  process_data() ;

  write_basis_functions() ;

  end = clock() ;

  printf("\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;
  fprintf(log_file, "\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  fclose(log_file) ;

#if USE_MPI == 1
  MPI_Finalize() ; 
#endif

  return 0; 
}

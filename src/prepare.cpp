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
 * If reaction types are known (know_reactions_flag = 1), then:
 *
 * reactions_in_channel:	indices of reactions in each reaction channel 
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

  // assign reactions to channels
  if (know_reactions_flag == 1)
  {
    reactions_in_channel.resize(channel_num) ;
    for (int i = 0; i < R; i ++)
    {
      if (channel_to_idx.count(vvec[i]) == 1)
      {
	// determine channel index of reaction i
	idx = channel_to_idx[vvec[i]] ;
	// add to reaction list of the channel 
	reactions_in_channel[idx].push_back(i) ;
      }       
      else // something wrong ...
      {
	printf("\nWarning: reaction %d doesn't correspond to any channel found in trajectories\n", i) ;
	fprintf(log_file, "\nWarning: reaction %d doesn't correspond to any channel found in trajectories\n", i) ;
      }
    }

    // make sure that each channel should contain at least one reaction
    for (int i = 0 ; i < channel_num; i ++)
      if ( reactions_in_channel[i].size() == 0 )
      {
	printf("Error: channel %d doesn't contain any reactions! Please check the trajectory data!\n", i) ;
	fprintf(log_file, "Error: channel %d doesn't contain any reactions! Please check the trajectory data!\n", i) ;
	exit(1);
      }
  }
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

  printf("\nInformation of channels are written into file: %s\n\n", buf) ;
  fprintf(log_file, "\nInformation of channels are written into file: %s\n\n", buf) ;

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

  printf("\nIn total:\n   %d reactions,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;
  fprintf(log_file, "\nIn total:\n   %d reactions,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;

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
 * (Task 1) If reaction types are known (know_reactions_flag=1), then the task of the
 * code ./sparse_learning is to learn the rate constants of reactions. 
 * In this case, each reaction corresponds to one basis function, and therefore 
 * there will be in total R basis functions. 
 *
 * The correspondence between reactions and basis functions are as follows:
 *
 *    Reaction    |   Basis
 * 1. A  -> *     |    x
 * 2. A+B -> *    |    x*y
 * 3. 2A -> *     |    x(x-1)
 *
 * In the above, A, B are two different species, and x, y are their copy-numbers.
 *
 * The (sparse) weights of rate constants are set to 0 in the output file.
 *
 * (Task 2) If reaction types are unknown (know_reactions_flag=0), the number of basis
 * functions (num_basis) depends on the value of poly_order. The task of 
 * the code ./sparse_learning is to approximate propensity functions of all
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

  /*
   * Record vectors of (non-repeated) reactants among the R reactions
   *
   * It will only be used in Task 1.
   *
   */
  map<vector<int>, int> vvec_in_map ;

  /* Records the index of basis function for each of R reactions.
   *
   * It will only be used in Task 1.
   */
  vector<int> basis_idx_of_reactions ;
  basis_idx_of_reactions.resize(R) ;

  // initialize the vector before including any basis functions
  basis_vec.resize(0) ;

  // if reaction types are known 
  if (know_reactions_flag == 1)
  {
    // (Task 1) in this case, we use the propensity functions of R reactions (excluding the repeated ones) as basis functions
    for (int i = 0; i < R; i ++)
      // if the same vector of reactants has not be recorded
      if ( vvec_in_map.count(vvec_in[i]) == 0 ) 
	{
	  reactant_idx.resize(0) ;
	  for (int j = 0; j < n; j ++)
	    for (int j1 = 0; j1 < vvec_in[i][j] ; j1 ++)
	      reactant_idx.push_back(j) ;
	    
	  // add this new basis 
	  basis_vec.push_back( reactant_idx ) ;

	  // record its index
	  vvec_in_map[ vvec_in[i] ] = basis_vec.size() - 1 ;
	  basis_idx_of_reactions[i] = basis_vec.size() - 1 ;
	}
      else // if the same basis has already occured
	basis_idx_of_reactions[i] = vvec_in_map[ vvec_in[i] ] ;

    num_basis = vvec_in_map.size() ; 
  }
  else // (Task 2)
  {
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
      num_basis = n * (3 + n) / 2 ;   
    } else 
    {
      num_basis = n ; 
    }
  }

  // initialization
  omega_vec.resize(channel_num) ;
  for (int i = 0 ; i < channel_num ; i ++)
    omega_vec[i].resize( num_basis, 0.0 ) ;

  // print basis functions into file
  printf("Number of basis functions = %d\n", num_basis) ;
  fprintf(log_file, "Number of basis functions = %d\n", num_basis) ;

  ofstream out_file ;
  char buf[100] ;
  sprintf( buf, "./output/basis_funct_info.txt" ) ;

  printf("Basis functions are written into the file: %s\n", buf) ;
  fprintf(log_file, "Basis functions are written into the file: %s\n", buf) ;
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

  out_file << channel_num << endl ;

  if (know_reactions_flag == 1) // (Task 1) if we know reaction structures
  {
    // output number of basis functions for each channel i, in this case, the
    // number of reactions in each channel 
    for ( int i = 0 ; i < channel_num ; i ++ )
      out_file << reactions_in_channel[i].size() << ' ' ;
    out_file << endl ;

    int ridx ;
    // reaction indices for each channel 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < reactions_in_channel[i].size() ; j ++)
      {
	// get the index of reaction first
	ridx = reactions_in_channel[i][j] ;
	// then output the index of the basis function that corresponds to this reaction
	out_file << basis_idx_of_reactions[ridx] << ' ' ;
      }
      out_file << endl ;
    }

    // the sparsity weight of each parameter equals to 0
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < reactions_in_channel[i].size() ; j ++)
	out_file << 0.0 << ' ' ;
      out_file << endl ;
    }

    // initial value for each basis functions 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < reactions_in_channel[i].size() ; j ++)
	out_file << 1.0 << ' ' ;
      out_file << endl ;
    }
  } else  // (Task 2) unknown reaction types
  {
    // output number of basis functions for each channel i, each channel used
    // all basis functions by default
    for ( int i = 0 ; i < channel_num ; i ++ )
      out_file << num_basis << ' ' ;
    out_file << endl ;

    // basis functions used for each channel 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < num_basis ; j ++)
	out_file << j << ' ' ;
      out_file << endl ;
    }

    // output weight of each parameter 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < num_basis ; j ++)
	out_file << 1.0 << ' ' ;
      out_file << endl ;
    }

    // initial value of parameters for each basis functions 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < num_basis ; j ++)
	out_file << 0.0 << ' ' ;
      out_file << endl ;
    }
  }

  out_file.close() ;
}

int main ( int argc, char * argv[] ) 
{
#if USE_MPI == 1
  MPI_Init(&argc, &argv) ;
#endif

  clock_t start , end ;
  char buf[30]; 

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

  output_channel_info() ;

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

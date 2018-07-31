#include "ssa.h"

// find all channels in the N_traj trajectories data
void find_channels_in_traj(vector<vector<vector<int> > > & traj_data)
{
  map<vector<int>, int> channel_to_idx ;
  set<vector<int> > channel_vec_set ;
  vector<int> vec_change ;
  int idx ;

  vec_change.resize(dim) ;
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
      for (int j = 0; j < dim; j ++)
	vec_change[j] = traj_data[traj_idx][i+1][j] - traj_data[traj_idx][i][j] ;

      // add a new channel to the set, if the change vector doesn't belong to the existing ones
      if ( channel_vec_set.count(vec_change) == 0 ) channel_vec_set.insert(vec_change) ;
    }
  }

  // channel index starts from 1
  // create list of channels 
  channel_num = 0 ;
  for (set<vector<int> >::iterator it = channel_vec_set.begin() ; it != channel_vec_set.end() ; it ++)
    {
      channel_num ++ ;
      channel_to_idx[*it] = channel_num ;
      channel_list.push_back(*it) ;
    }

  for (int traj_idx = 0 ; traj_idx < N_traj ; traj_idx ++)
  {
    Mi_in_traj[traj_idx].resize(channel_num, 0) ;
    // for each reaction in the trajectory 
    for (int i = 0; i < num_state_in_traj[traj_idx]-1; i ++)
    {
      // compute change vector
      for (int j = 0; j < dim; j ++)
	vec_change[j] = traj_data[traj_idx][i+1][j] - traj_data[traj_idx][i][j] ;

      // determine the index for each reaction
      idx = channel_to_idx[vec_change] ;
      channel_idx_in_traj[traj_idx][i] = idx ;
      // increase counter 
      Mi_in_traj[traj_idx][idx - 1] ++ ;
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
	reactions_in_channel[idx-1].push_back(i) ;
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

// output the information of channels to file 
void output_channel_info()
{
  ofstream out_file ;
  char buf[100] ;
  sprintf( buf, "./output/channel_info.txt") ;

  printf("\nInformation of channels are stored to file: %s\n", buf) ;
  fprintf(log_file, "\nInformation of channels are stored to file: %s\n", buf) ;

  out_file.open(buf) ;

  // make sure the file is open successfully
  if ( ! out_file.is_open() )
    {
      printf("Error: can not open file : %s \n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s \n\n", buf) ;
      exit(1) ;
    }

  out_file << channel_num << ' ' << dim << endl ;

  for (int i = 0; i < channel_num; i ++)
  {
    for (int j = 0; j < dim ; j++)
      out_file << channel_list[i][j] << ' ';
    out_file << endl ;
  }

  out_file << endl ;

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

  printf("No. of trajectories = %d\n\n", N_traj) ;
  fprintf(log_file , "No. of trajectories = %d\n\n", N_traj) ;

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

    in_file >> dim ;
    c_state.resize(dim) ;

    num_state_in_traj[traj_idx] = 0 ;

    // read the file line by line
    while ( getline(in_file, line) )
    {
      if (line.find_first_not_of(' ') != string::npos) 
      {
	istringstream iss(line) ;
	iss >> t_now ;
	for (int i = 0 ; i < dim ; i ++)
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

  printf("\n%d reactions in total,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;
  fprintf(log_file, "\n%d reactions in total,\t%d reaction channels\n\n", tot_num_reaction_in_traj, channel_num) ;

  for (int i =0 ; i < channel_num ; i ++)
  {
    printf("Change vector of the %dth channel :   ", i+1) ; 
    fprintf(log_file, "Change vector of the %dth channel :   ", i+1) ; 
    for (int j = 0 ; j < dim ; j ++)
    {
      cout << channel_list[i][j] << ' ' ;
      fprintf( log_file, "%d ", channel_list[i][j] );
    }
    cout << endl ;
    fprintf(log_file, "\n") ;
  }

  printf("Occurrence of each reaction channels in these %d trajectories : [%d", N_traj, Mi_in_all_traj[0]) ;
  fprintf(log_file, "Occurrence of each reaction channels in these %d trajectories : [%d", N_traj, Mi_in_all_traj[0]) ;
  for (int i = 1 ; i < channel_num ; i ++)
  {
    printf(", %d", Mi_in_all_traj[i]) ;
    fprintf(log_file, ", %d", Mi_in_all_traj[i]) ;
  }
  printf("]\n") ;
  fprintf(log_file, "]\n") ;
}

// write basis functions to file
void write_basis_functions() 
{
  vector<int> reactant_idx ;

  reactant_idx.resize(1) ; 
  // if reaction types are known, then only the corresponding propensity
  // functions as basis functions
  if (know_reactions_flag == 1)
  {
    // in this case, we use the R propensity functions basis functions
    for (int i = 0; i < R; i ++)
    {
      reactant_idx.resize(0) ;
      for (int j = 0; j < dim ; j ++)
	for (int j1 = 0; j1 < vvec_in[i][j] ; j1 ++)
	  reactant_idx.push_back(j) ;
	
      basis_vec.push_back( reactant_idx ) ;
    }
    num_basis = R ; 
  }
  else 
  {
    // include polynomial functions of order 1 (linear) as basis functions
    for (int i =0 ; i < dim ; i++)
    {
      reactant_idx[0] = i ;
      basis_vec.push_back( reactant_idx ) ;
    }

    // include quadratic functions as basis functions
    if (poly_order == 2)
    {
      reactant_idx.resize(2) ; 
      for (int i=0 ; i < dim ; i++)
	for (int j=i ; j < dim ; j++)
	{
	  reactant_idx[0] = i ;
	  reactant_idx[1] = j ;
	  basis_vec.push_back( reactant_idx ) ;
	}
      num_basis = dim * (3 + dim) / 2 ;   
    } else 
    {
      num_basis = dim ; 
    }
  }

  // initialization
  omega_vec.resize(channel_num) ;
  for (int i = 0 ; i < channel_num ; i ++)
    omega_vec[i].resize( num_basis, 0.0 ) ;

  // print basis functions into file
  printf("No. of basis functions = %d\n", num_basis) ;
  fprintf(log_file, "No. of basis functions = %d\n", num_basis) ;

  ofstream out_file ;
  char buf[100] ;
  sprintf( buf, "./output/basis_funct_info.txt" ) ;

  printf("Basis functions are stored into the file: %s\n", buf) ;
  fprintf(log_file, "Basis functions are stored into the file: %s\n", buf) ;
  out_file.open(buf) ;
  if ( ! out_file.is_open() )
    {
      printf("Error: can not open output file : %s. \n\n", buf) ;
      fprintf(log_file, "Error: can not open output file : %s. \n\n", buf) ;
      exit(1) ;
    }

  out_file << dim << ' ' << num_basis << endl ;
  for (int i = 0 ; i < num_basis ; i ++)
  {
    out_file << basis_vec[i].size() << ' ' ;
    for (int j = 0 ; j < basis_vec[i].size() ; j ++)
      out_file << basis_vec[i][j] << ' ' ;
    out_file << endl ;
  }
  out_file << endl ;

  out_file << channel_num << endl ;

  if (know_reactions_flag == 1) // if we know reaction structures, this is the Task 1
  {
    // output number of basis functions for each channel i, in this case, the
    // number of reactions in each channel 
    for ( int i = 0 ; i < channel_num ; i ++ )
      out_file << reactions_in_channel[i].size() << ' ' ;
    out_file << endl ;

    // reaction indices for each channel 
    for ( int i = 0 ; i < channel_num ; i ++ )
    {
      for (int j = 0 ; j < reactions_in_channel[i].size() ; j ++)
	out_file << reactions_in_channel[i][j] << ' ' ;
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
  } else  // this is the Task 2 (unkown reaction types)
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
  MPI_Init(&argc, &argv) ;

  clock_t start , end ;
  char buf[30]; 

  sprintf(buf, "./log/prepare.log") ;
  if (init(buf) < 0) return -1 ;

  // this program should be ran sequentially
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

  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;
  fprintf(log_file, "\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  fclose(log_file) ;

  MPI_Finalize() ; 

  return 0; 
}

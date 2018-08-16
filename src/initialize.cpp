#include "sparse_learning.h"

int read_config() ;

/* 
 * Initialize the seed of random numbers
 * The seeds depend on both the time of the machine and the index of the local processor 
 *
 */
void init_rand_generator()
{
  long is1, is2 ;

  char phrase[100] ;
  sprintf( phrase, "%ld", time(NULL) + mpi_rank ) ;

  phrtsd(phrase, &is1, &is2) ;
  setall(is1, is2) ;

  // test 10 random numbers and print to log file 
  if (mpi_rank == 0)
  {
    fprintf(log_file, "\nGenerate 10 random numbers (uniform on [0,1]) :\n ");

    for (int i = 0 ; i < 10 ; i ++)
      fprintf( log_file, "%.8f\t", ranf() ) ;

    fprintf(log_file, "\n");
  }
}

/*
 * Read the information of reaction network from file
 */

void read_reactions()
{
  char buf[100] ;

  sprintf( buf, "./input/cr.txt") ;
  ifstream cr( buf ) ; 

  if (mpi_rank == 0)
  {
    printf("\nReading reaction information from the file: %s\n", buf) ;
    fprintf(log_file, "\nReading reaction information from the file: %s\n", buf) ;
  }

  // check whether the file is successfully open 
  if ( (mpi_rank == 0) && (cr.is_open() == 0) ) 
    {
      printf("Error: can not open file : %s.\n\n", buf) ;
      fprintf(log_file, "Error: can not open file : %s.\n\n", buf) ;
      exit(1) ;
    }

  // read reaction number(R) and state dimension (n) 
  cr >> R >> n ;

  vvec_in.resize(R) ; vvec_out.resize(R) ; vvec.resize(R) ;
  init_state.resize(n) ;
  v_neg_idx_vec.resize(R) ;
  reactant_num.resize(R) ;

  // read the initial state
  for (int j = 0 ; j < n; j ++)
    cr >> init_state[j] ;

  // read the info of each reaction 
  for (int i = 0 ; i < R ; i ++)
  {
    vvec_in[i].resize(n) ; vvec_out[i].resize(n) ; vvec[i].resize(n) ;
    reactant_num[i] = 0 ;
    // read in vector of each reaction 
    for (int j = 0 ; j < n ; j ++)
    {
      cr >> vvec_in[i][j] ;
      reactant_num[i] += vvec_in[i][j] ;
    }
    // read out vector of each reaction 
    for (int j = 0 ; j < n ; j ++)
      cr >> vvec_out[i][j] ;

    // compute change vector of each reaction 
    for (int j = 0 ; j < n ; j ++)
    {
      vvec[i][j] = vvec_out[i][j] - vvec_in[i][j] ;
      if (vvec[i][j] < 0) v_neg_idx_vec[i].push_back(j) ;
    }
  }

  // load rate constant kappa[i] for ith reaction
  kappa_vec.resize(R) ;
  for (int i = 0 ; i < R ; i ++)
    cr >> kappa_vec[i] ;

  cr.close() ;
}

int init(char * log_file_name )
{

#if USE_MPI == 1
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank ) ;
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size ) ;
#else
  mpi_rank = 0 ;
  mpi_size = 1 ;
#endif

  if (mpi_rank == 0)
  {
    char buf[30] ;
    // create directory for log file 
    sprintf( buf, "./log/") ;
    if (dir_check(buf) != 0) 
    {
      printf ("\nError: directory ./log/ can not be created. Fix the error manually!\n\n") ;
      exit(1) ;
    }

    // prepare log file 
    log_file = fopen(log_file_name, "w") ;
    if ( log_file == NULL )
      {
	printf( "Error: can not open log file : %s. \n\n", log_file_name ) ;
	exit(1) ;
      }
  }

  init_rand_generator() ;

  // read parameters 
  if (read_config() < 0) return -1 ;

  // distribute N_traj trajectories among processors
  local_N_traj = N_traj / mpi_size ; 
  local_traj_start_idx = (N_traj / mpi_size) * mpi_rank ;
  if ( mpi_rank < N_traj % mpi_size ) 
  {
    local_N_traj ++ ;
    local_traj_start_idx += mpi_rank ;
  }
  else local_traj_start_idx += N_traj % mpi_size ;

  // if reaction types are known, then read reactions from file 
  if (know_reactions_flag == 1) read_reactions() ;

  if (mpi_rank == 0)
  { // print information 
    if (know_reactions_flag == 1)
    {
      printf( "\nNumber of reactions = %d,\tDim = %d\n", R, n) ;
      fprintf( log_file, "\nNumber of reactions = %d,\tDim = %d\n", R, n) ;
    }

    printf( "\nNumber of processors = %d\n\n", mpi_size ) ;
    fprintf( log_file, "\nNumber of processors = %d\n\n", mpi_size ) ;
  }


  return 0 ;
}

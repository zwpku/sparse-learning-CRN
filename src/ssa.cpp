#include "sparse_learning.h"

/*
 * Generate trajectories using SSA method.
 *
 * Each local processor generates local_N_traj trajectories.
 *
 * In total, N_traj trajectories of length (time) T will be saved under
 * directory ./traj_data/
 *
 */

void gen_data() 
{
  char buf[100] ;
  ofstream out_file ;
  double t_now , tau ;
  vector<int> c_state , next_state ;

  next_state.resize(n) ;

  // each proccessor generates several (=local_N_traj) trajectories
  for ( int j = 0 ; j < local_N_traj ; j ++ )
    {
      // open file 
      sprintf( buf, "./traj_data/traj_%d.txt", j + local_traj_start_idx ) ;
      out_file.open(buf) ;

      // make sure no error occurs 
      if ( ! out_file.is_open() )
      {
	printf("Error: can not open file : %s. \n\n", buf) ;
	fprintf(log_file, "Error: can not open file : %s. \n\n", buf) ;
	exit(1) ;
      }

      out_file << n << endl ;

      t_now = 0 ; 
      c_state = init_state ;

      while ( t_now < T )  // up to time T 
      {
	// generate next state (next_state) from current state (c_state)
	tau = ssa(t_now, c_state, next_state) ;
	// no more reactions 
	if (tau + t_now > T) tau = T - t_now ;

	out_file << t_now << ' ' ;
	for (int i = 0 ; i < n ; i ++)
	  out_file << c_state[i] << ' ' ;
	out_file << tau << endl ;

	// update state and time
	c_state = next_state ; 
	t_now += tau ; 
      }

      out_file.close() ;
    }
}

int main ( int argc, char * argv[] ) 
{
  MPI_Init(&argc, &argv) ;

  clock_t start , end ;
  char buf[30]; 

  // prepare the file for printing log information
  sprintf(buf, "./log/ssa.log") ;
  if (init(buf) < 0) return -1 ;

  // for SSA, reactions information is need, therefore the input parameter is ignored
  if (know_reactions_flag == 0) 
  {
    know_reactions_flag = 1;

    if (mpi_rank == 0)
    {
      printf("Variable know_reaction_flag is reset to 1\n");
      fprintf(log_file, "Variable know_reaction_flag is reset to 1\n");
    }
  }

  if (mpi_rank == 0)
  {
    // make sure the directory is ready
    sprintf(buf, "./traj_data") ;
    if (dir_check(buf) != 0) 
    {
      printf ("\nError: directory %s can not be created. Fix the error manually!\n\n", buf) ;
      fprintf (log_file, "\nError: directory %s can not be created. Fix the error manually!\n\n", buf) ;
      exit(1) ;
    }
  }

  start = clock() ;

  gen_data() ;

  end = clock() ;

  if (mpi_rank == 0)
  {
    printf(" %d trajectories have been generated under directory: ./traj_data/\n", N_traj) ;
    printf(" Runtime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

    fprintf(log_file, " %d trajectories have been generated under directory: ./traj_data/\n", N_traj) ;
    fprintf(log_file, " Runtime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

    fclose(log_file) ;
  }

  MPI_Finalize() ; 

  return 0; 
}

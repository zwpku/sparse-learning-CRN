#include "sparse_learning.h"
#include <iostream>
#include "libconfig.h++"
 
using namespace libconfig;

template <class T>
int read_value(Config & cfg, string  name, T & val)
{
  try
  {
    cfg.lookupValue(name, val);
  }
  catch(const SettingNotFoundException &nfex)
  {
    if (mpi_rank == 0)
      cerr << "No " << name << " setting in configuration file." << endl;
    return -1;
  }
  return 0;
}

/*
 *
 * Read parameters from the file: sparse_learning.cfg 
 *
 */

int read_config() 
{
  Config cfg;

  try {
    cfg.readFile("sparse_learning.cfg");
  }
  catch (const FileIOException &fioex)
  {
    if (mpi_rank == 0)
      cerr << "I/O error while reading file." << endl;
    return -1;
  }
  catch(const ParseException &pex)
  {
    if (mpi_rank == 0)
      cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << endl ;
    return -1;
  }

  if ( read_value(cfg, string("T"), T ) < 0 )
    return -1;

  if ( read_value(cfg, string("N_traj"), N_traj) < 0 )
    return -1;

  if ( read_value(cfg, string("poly_order"), poly_order) < 0 )
    return -1;

  if ( read_value(cfg, string("xx_basis_flag"), xx_basis_flag) < 0 )
    return -1;

  if ( read_value(cfg, string("regular_lambda"), regular_lambda) < 0 )
    return -1;

  if ( read_value(cfg, string("epsL1_flag"), epsL1_flag) < 0 )
    return -1;

  if ( read_value(cfg, string("solver_id"), solver_id) < 0 )
    return -1;

  if ( read_value(cfg, string("grad_dt"), grad_dt) < 0 )
    return -1;

  if ( read_value(cfg, string("eps"), eps) < 0 )
    return -1;

  if ( read_value(cfg, string("cost_stop_tol"), cost_stop_tol) < 0 )
    return -1;

  if ( read_value(cfg, string("tot_step"), tot_step) < 0 )
    return -1;

  if ( read_value(cfg, string("max_step_since_prev_min_cost"), max_step_since_prev_min_cost) < 0 )
    return -1;

  if ( read_value(cfg, string("num_record_tail_cost"), num_record_tail_cost) < 0 )
    return -1;

  if ( read_value(cfg, string("delta"), delta) < 0 )
    return -1;

  if ( read_value(cfg, string("g_cut"), g_cut) < 0 )
    return -1;

  if ( read_value(cfg, string("output_interval"), output_interval) < 0 )
    return -1;

  if ( read_value(cfg, string("know_reactions_flag"), know_reactions_flag) < 0 )
    return -1;

  /*
   * If reaction types are known, then use basis x*(x-1) instead of x*x.
   *
   * See the document in the file prepare.cpp for details.
   */
  if (know_reactions_flag == 1) xx_basis_flag = 0 ;

  return 0 ;
}


#include "pgpr_ppic.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[])
{
  MPI_Init( NULL, NULL );
  pgpr_parse cfg(CFGDEMO, argc, argv);
  Int mnum; //the number of machines
  MPI_Comm_size( MPI_COMM_WORLD, &mnum );
  if( mnum != cfg.param_demo.blocks ) {
    MPI_Finalize();
    throw("There are not enough machines available.\n");
  }

  pgpr_ppic gp(cfg.param_demo.hyperf);
  gp.regress(cfg.param_demo.trainf, cfg.param_demo.testf,
             cfg.param_demo.supportf, cfg.param_demo.blocks);

  Int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if( rank == 0 ) {
    pmsg(LEV_PRG, stdout, "pPIC  |");
    gp.outputRst(cfg.param_demo.outf);
  }
  MPI_Finalize();

}

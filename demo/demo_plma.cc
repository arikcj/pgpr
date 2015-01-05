#include "pgpr_plma.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[]){

  MPI_Init(NULL, NULL);
  pgpr_parse cfg(CFGDEMO, argc, argv);

  Int mnum; //the number of machines
  MPI_Comm_size(MPI_COMM_WORLD, &mnum);
  
  if (mnum != cfg.param_demo.blocks ) {
    MPI_Finalize();
    throw("There are not enough machines availabe.\n");
  }


  pgpr_plma gp(cfg.param_demo.hyperf, cfg.param_demo.trainf, cfg.param_demo.testf,
			   cfg.param_demo.supportf, cfg.param_demo.bandwidth, cfg.param_demo.blocks);
  gp.regress();

  Int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0){
    pmsg(LEV_PRG, stdout, "pLMA  |");
    gp.outputRst(cfg.param_demo.outf);
  }

  MPI_Finalize();
}

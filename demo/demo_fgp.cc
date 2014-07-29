#include "pgpr_fgp.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[])
{
  pgpr_parse cfg(CFGDEMO, argc, argv);
  pgpr_fgp gp(cfg.param_demo.hyperf);
  pmsg(LEV_PRG, stdout, " FGP  |");
  gp.regress(cfg.param_demo.trainf, cfg.param_demo.testf);
  gp.outputRst(cfg.param_demo.outf);

}

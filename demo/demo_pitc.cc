#include "pgpr_pitc.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[])
{
  pgpr_parse cfg(CFGDEMO, argc, argv);
  pgpr_pitc gp(cfg.param_demo.hyperf);
  pmsg(LEV_PRG, stdout, " PITC |");
  gp.regress(cfg.param_demo.trainf, cfg.param_demo.testf,
             cfg.param_demo.supportf, cfg.param_demo.blocks);
  gp.outputRst(cfg.param_demo.outf);

}

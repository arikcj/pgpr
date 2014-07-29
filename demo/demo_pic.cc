#include "pgpr_pic.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[])
{
  pgpr_parse cfg(CFGDEMO, argc, argv);
  pgpr_pic gp(cfg.param_demo.hyperf);
  pmsg(LEV_PRG, stdout, " PIC  |");
  gp.regress(cfg.param_demo.trainf, cfg.param_demo.testf,
             cfg.param_demo.supportf,cfg.param_demo.blocks);
  gp.outputRst(cfg.param_demo.outf);

}

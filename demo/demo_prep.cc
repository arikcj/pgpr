#include "pgpr_util.h"
#include "pgpr_parse.h"
#include "pgpr_data.h"

int main(int argc, char * argv[])
{
  pgpr_parse cfg(CFGPREP, argc, argv);
  Int mach_num = cfg.param_prep.mach_num;
  Int blk_size = cfg.param_prep.blk_size;
  Int ddim = cfg.m_xyz.ncols();
  Int dsize = cfg.m_xyz.nrows();

  saveHyper(cfg.param_prep.hyperf, cfg.v_hyp, cfg.in_dim);

  pgpr_data data(cfg, mach_num, blk_size);
  //Testing set
  Int percent = cfg.param_prep.percent;
  Doub test_size = floor(dsize / percent);
  pmsg(LEV_DBG, stdout, "Intialize test set: %d \% [%d] of the domain\n", percent, test_size);
  Mdoub test_set;
  data.getTestSet(test_size, test_set);

  //save testing set
  saveData(cfg.param_prep.testf, test_set);

  //Training set
  pmsg(LEV_DBG, stdout, "Initialize prior data\n");
  data.getTrainSet(ALGO1);

  Int ds = 0;
  Mdoub train_data(mach_num * blk_size, ddim);
  for(Int i = 0; i < mach_num; i++) {
    for (Int j = 0; j < blk_size; j++) {
      Int pos = (Int) data.m_sample[j][i];
      takeSamp(cfg.m_xyz[pos], train_data[ds], ddim);
      ds++;
    }
  }
  saveData(cfg.param_prep.trainf, train_data);

//Support set
  if (cfg.param_prep.mode == 1) {
    Int support = cfg.param_prep.support;
    Mdoub support_set(support, ddim);
    pmsg(LEV_DBG, stdout, "Initialize support set: %d\n", support);
    data.selMaxVar(cfg.m_xyz, dsize, support, support_set);

    saveData(cfg.param_prep.supportf, support_set);

  }

}

/** @file pgpr_parse.h
 * 
 * @brief This file provides functionalities to parse commandlines of
 * different applications and domain files.
 * 
 * @detail In this version, this class support parsing two type of
 * applications: 1. preparing data for experiment; 2. demonstration of
 * regression algorithms. In addition, this class can load and parse
 * domain files which includs a file storing the domain inputs and outputs 
 * (features and targets), and a configuration file containing basic information 
 * of the domain.
 * 
 * @author CHEN jie, arik.cj@gmail.com
 * 
 * @version 1.0
 */

#ifndef _PGPR_PARSE_H_
#define _PGPR_PARSE_H_
#include "pgpr_type.h"
#include "pgpr_util.h"
#define DOM 0  // Domain size
#define INP 1  // Input dimension
#define OUT 2  // Output dimension
#define MNE 3  // Maximum size of neighbor nodes
#define HYP 4  // Number of hyper
#define PARANUM 5 // No. of parameters in first line of cfg file
#define NAMELEN 128 // Length of a file name
#define CFGDEMO 1 // Demonstration application 
#define CFGPREP 2 // Experimental data preparation application
/**
 * @struct t_command_prep	
 * @brief Information parsed from commandline of application that
 * prepares the experimental data.
 **/
struct t_command_prep {
	/// Running mode of program
  Int mode;			
	/// Random seed
  Int seed;				
	/// No. of blocks/machines
  Int mach_num;				
	/// Data size in each block/machine
  Int blk_size;				
	/// File containing the training data
  Char trainf[NAMELEN];
	/// Percentage of domain inputs used for testing
  Int percent;			
	/// File containing the test set
  Char testf[NAMELEN];
	/// File containing the hyperparameters
  Char hyperf[NAMELEN];
	/// Size of support set
  Int support;				
	/// File containing the support set (used by PITC/pPITC/PIC/pPIC)
  Char supportf[NAMELEN];
};

/**	
 * @struct t_command_demo	
 * @brief Information parsed from commandline of application that
 * demonstrates the regression algorithms. 
 **/
struct t_command_demo {
	/// Running mode (reserved)
  Int mode;								
	/// Number of data blocks (used by PITC/pPITC/PIC/pPIC)
  Int blocks;
	/// Size of the reduced rank (used by PICF-based GP)
  Int rank;							
	/// File containing the results
  Char outf[NAMELEN];		
	/// File containing the training data
  Char trainf[NAMELEN];	
	/// File containing the test set
  Char testf[NAMELEN];
	/// File containing the hyperparameters
  Char hyperf[NAMELEN];	
	/// File containing the support set (used by PITC/pPITC/PIC/pPIC)
  Char supportf[NAMELEN];	
};

/**
 * @class pgpr_parse
 * @brief This class parses domain data files, configuration file and commandline of
 * different applications
 */
class pgpr_parse
{
public:
	/// Vector of parameters loaded from cfg file
  Vdoub v_param;
	/// Vector of hyperparameters loaded from cfg file
  Vdoub v_hyp;	
	/// @brief Matrix storing domain data where no. of rows indicates the size of
	/// domain and no. of columns is the dimension of inputs and outputs
  Mdoub m_xyz;	
	/// Dimension of inputs
  Int in_dim;
	/// Dimension of outputs
  Int out_dim;
	/// @brief Connectivity of the domain data: each row stores the indices of a set of 
	///data points that connects to the corresponding data points.
  Mint  m_connect;
	/// Information extracted from the commandline of a demo application
  t_command_demo param_demo;
	/// @brief Information extracted from the commandline of experimental data
	/// preparing application
  t_command_prep param_prep;

	/**
	 *	@brief Construction function
	 *	@param[in] md Different applications: CFGDEMO - demo application; CFGPREP - data preparing
	 *	application
	 *	@param[in] argc count of arguments 
	 *	@param[in] argv an array of pointers to the list of arguments
	 **/
  pgpr_parse(Int md, Int argc, Char * argv[]) {
    switch(md) {
    case CFGDEMO:
      parse_commandline_demo(argc, argv);
      break;
    case CFGPREP:
      parse_commandline_prep(argc, argv);
      break;
    default:
      throw("The commandline is not supported\n");
    }
  }

private:

	/**
	 *	@brief Parse a domain data file 
	 *	@param[in] file Name of a domain data file 
	 **/
  void parse_domfile(Char * file) {
    FILE * fp;
    Int ret;
    fp = fopen(file, "r");
    if(fp == NULL) {
      throw("Fail to open file\n");
    }
    Int d_xz = (Int)v_param[INP] + (Int)v_param[OUT];
    m_xyz.resize((Int)v_param[DOM], d_xz);
    m_connect.resize((Int)v_param[DOM], (Int)v_param[MNE] + 1);
    for(Int i = 0; i < v_param[DOM]; i++) {
      for(Int j = 0; j < d_xz; j++) {
        Doub tmp;
        ret = fscanf(fp, "%lf", &tmp);
        m_xyz[i][j] = tmp;
      }
      Int num_ne;
      ret = fscanf(fp, "%d", &num_ne);
      m_connect[i][0] = num_ne;
      for(Int j = 1; j <= num_ne; j++) {
        Int tmp;
        ret = fscanf(fp, "%d", &tmp);
        m_connect[i][j] = tmp;
      }
      ret = fscanf(fp, "\n");

    }
    fclose(fp);
  }
  
	/**
	 *	@brief Parse a configuration file 
	 *	@param[in] file Name of a configuration file corresponding to a
	 *	domain 
	 **/
  void parse_cfgfile(Char * file) {
    FILE * fp;
    Int ret;
    fp = fopen(file, "r");
    if(fp == NULL) {
      throw("Fail to open file\n");
    }
    //cofig of the domain data (#parameter is predefined)
    ret = fscanf(fp, "dom:");
    for(Int i = 0; i < PARANUM; i++) {
      Doub tmp;
      ret = fscanf(fp, "%lf", &tmp);
      v_param[i] = tmp;
      pmsg(LEV_DBG, stdout, "%lf ", v_param[i]);
    }
    ret = fscanf(fp, "\n");
    pmsg(LEV_DBG, stdout, "\n");
    //config of hyper-parameters
    ret = fscanf(fp, "hyp: ");
    v_hyp.resize((Int)v_param[HYP]);
    for(Int i = 0; i < v_hyp.size(); i++) {
      Doub tmp;
      ret = fscanf(fp, "%lf", &tmp);
      pmsg(LEV_DBG, stdout, "%lf ", tmp);
      v_hyp[i] = tmp;
    }
    ret = fscanf(fp, "\n");
    pmsg(LEV_DBG, stdout, "\n");
    in_dim = (Int)v_param[INP];
    out_dim = (Int)v_param[OUT];
    fclose(fp);
  }
	
	/**
	 *	@brief Parse the commandline of a data preparing application
	 *	@param[in] argc count of arguments 
	 *	@param[in] argv an array of pointers to the list of arguments
	 **/
  Int parse_commandline_prep(Int argc, Char * argv[]) {
    v_param.resize(PARANUM);
    if (argc < 8) {
      pmsg(LEV_PRG, stdout, "Usage: <app> -dataset <domain> -mode <mode>\
-seed <seed> -machine <no. of blocks/machines> -blocksize <size of \
each block> -support <size of support set> -percent <percentage of \
the domain data used as test set> \n");
			throw("Number of arguments in commandline is not enough.\n");
    } else { // if we got enough parameters...
      Char * data_key;
      Char dom_file[NAMELEN];
      Char cfg_file[NAMELEN];
      Char * out_key;
      for (int i = 1; i < argc; i++) {
		  if (i + 1 != argc){
			  if (!strcmp(argv[i] , "-dataset")) {
				  data_key = argv[++i];
				  sprintf(dom_file, "%s.%s", data_key, "dom");
				  sprintf(cfg_file, "%s.%s", data_key, "cfg");

				  pmsg(LEV_DBG, stdout, "config:%s\n", cfg_file);
				  pmsg(LEV_DBG, stdout, "domain:%s\n", dom_file);
			  }else if (!strcmp(argv[i] , "-seed")) {
				  param_prep.seed = atoi(argv[++i]);
				  SRAND(param_prep.seed);
				  pmsg(LEV_DBG, stdout, "random seed:%d\n", param_prep.seed);
			  } else if (!strcmp(argv[i] , "-output")) {
				  out_key = argv[++i];
			  } else if (!strcmp(argv[i] , "-mode")) {
				  param_prep.mode = atoi(argv[++i]);
			  } else if (!strcmp(argv[i] , "-machine")) {
				  param_prep.mach_num = atoi(argv[++i]);
			  } else if (!strcmp(argv[i] , "-blocksize")) {
				  param_prep.blk_size = atoi(argv[++i]);
			  } else if (!strcmp(argv[i] , "-support")) {
				  param_prep.support = atoi(argv[++i]);
			  } else if (!strcmp(argv[i] , "-percent")) {
				  param_prep.percent = atoi(argv[++i]);
			  } else {
				  throw("Invalid arguments.\n");
			  }
		  }
	  }
      parse_cfgfile(cfg_file);
      parse_domfile(dom_file);
      sprintf(param_prep.trainf, "%s.%d.%s", out_key, param_prep.seed, "trn");
      sprintf(param_prep.testf, "%s.%d.%s", out_key, param_prep.seed, "tst");
      sprintf(param_prep.supportf, "%s.%d.%s", out_key, param_prep.seed, "spt");
      sprintf(param_prep.hyperf, "%s.%s", out_key, "hyp");
      pmsg(LEV_DBG, stdout, "training data:%s\n", param_prep.trainf);
      pmsg(LEV_DBG, stdout, "test data:%s\n", param_prep.testf);
      pmsg(LEV_DBG, stdout, "support set:%s\n", param_prep.supportf);
      pmsg(LEV_DBG, stdout, "hyperparameter:%s\n", param_prep.hyperf);
    }
    return SUCC;
  }
 
	/**
	 *	@brief Parse the commandline of a demonstration application
	 *	@param[in] argc count of arguments 
	 *	@param[in] argv an array of pointers to the list of arguments
	 *  @return SUCC
	 **/
  Int parse_commandline_demo(Int argc, Char * argv[]) {
    v_param.resize(PARANUM);
    if (argc < 2) {
      pmsg(LEV_PRG, stdout, "Usage: <app> -hyper <hyperparameters file> -in\
<an experimental case> -out <output file> -blocks <no. of blocks> -rank\
<reduced rank> \n");
			throw("Number of arguments in commandline is not enough.\n");
    } else { // if we got enough parameters...
      Char * data_key;
      Char dom_file[NAMELEN];
      Char cfg_file[NAMELEN];
      Char * in_key;
      for (int i = 1; i < argc; i++) {
		  if (i + 1 != argc){
			  if (!strcmp(argv[i] , "-hyper")) {
				  data_key = argv[++i];
				  sprintf(param_demo.hyperf, "%s", data_key);
			  } else if (!strcmp(argv[i] , "-in")) {
				  in_key = argv[++i];
				  sprintf(param_demo.trainf, "%s.trn", in_key);
				  sprintf(param_demo.testf, "%s.tst", in_key);
				  sprintf(param_demo.supportf, "%s.spt", in_key);
				  pmsg(LEV_DBG, stdout, "input:%s\n", param_demo.trainf);
				  pmsg(LEV_DBG, stdout, "input:%s\n", param_demo.testf);
				  pmsg(LEV_DBG, stdout, "input:%s\n", param_demo.supportf);
			  } else if (!strcmp(argv[i] , "-out")) {
				  Char *out_key = argv[++i];
				  sprintf(param_demo.outf, "%s.rst", out_key);
			  } else if (!strcmp(argv[i] , "-blocks")) {
				  param_demo.blocks = atoi(argv[++i]);
			  } else if (!strcmp(argv[i] , "-rank")) {
				  param_demo.rank = atoi(argv[++i]);
			  } else {
				  throw("Invalid arguments.\n");
			  }
		  }
	  }
    }
    return SUCC;
  }


};
#endif

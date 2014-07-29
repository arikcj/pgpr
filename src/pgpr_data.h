/** @file pgpr_data.h
 *  @brief This file provides class for maintaining and operating on a
 *  domain dataset, and class for training set, test set and support set
 *  selection from domain.
 *
 *  @author CHEN jie, arik.cj@gmail.com
 * 
 *  @version 1.0
 */

#ifndef _PGPR_DATA_H_
#define _PGPR_DATA_H_
#include "pgpr_type.h"
#include "pgpr_util.h"
#include "pgpr_parse.h"
#include "pgpr_fgp.h"
#define ALGO1 1 
#define ALGO2 2
//The index of an input starts from 1 
//The array index starts from 0
#define s2i(s) (s-1) 
#define i2s(i) (i+1) 
/**
 * @brief Status (e.g., modified or not) of a set of elements 
 **/
struct t_state {
  /// Status of a set of elements 
	Vbool v_st;
	/// Size of elements whose statuses are changed	
  Int ss;			

	t_state() {
    ss = 0;
  }
	
	/**
	 * @brief Constructor
	 * @param[in] ds size of the set
	 **/
  t_state(Int ds) {
    v_st.assign(ds, false);
    ss = 0;
  }

	/**
	 * @brief Initialize the structure
	 * @param[in] ds size of the set
	 **/
  void init(Int ds) {
    v_st.assign(ds, false);
  }

	/**
	 * @brief Get the status of a specific element
	 * @param[in] s Index of an element
	 * @return Binary status
	 **/
  inline Bool getSt(Int s) {
    return v_st[s2i(s)];
  }

	/**
	 * @brief Set the status of a specific element
	 * @param[in] s Index of an element
	 **/
  inline void upSt(Int s) {
    v_st[s2i(s)] = true;
    ss++;
  }

};
/** @class pgpr_domain
 *  @brief This class stores the domain data in a matrix manner
 */
class pgpr_domain
{
private:
  Mdoub attr; //attributs
  Mint ne;    //topology information
  Int deg;    //max degree 

public:
	///dimension of the domain (dimension of inputs and outputs).
  Int dd;     
	///size of the domain
  Int ds;
	
	///@brief Constructor
	///@param[in] n Matrix representing the topology of domain inputs
	///@param[in] a Matrix containing all domain data, both inputs and
	///outputs
  pgpr_domain(Mint n, Mdoub a): ne(n), attr(a) {
    dd = attr.ncols();
    ds = attr.nrows();
  }
	
	///@brief Compute the Euclidean distance between a pair of inputs
	///@param[in] ci index of an input
	///@param[in] cj index of an input
	///@return Value of Euclidean distance
  Doub getDist(Int ci, Int cj) {
    Doub sq_dist = 0;
    for(Int i = 0; i < dd - 1; i++) {
      sq_dist += SQR(attr[ci][i] - attr[cj][i]);
    }
    return SQRT(sq_dist);
  }

	///@brief Select randomly an index of a set, which is not selected
	///@param[in] st the status indicating if an element is selected or not
	///@return the index of an input 
	Int selRand(t_state *st) {
    Int s = -1;
    while(st->ss < ds) {
      s = RANI(ds - 1) + 1;
      if((st->getSt(s)) == 0 && (attr[s - 1][dd - 1] >= 0)) {
        st->upSt(s);
        break;
      }
    }
    return s;
  }

	///@brief Get the attributes of an input and its corresponding outputs
	///@param[in] s index of an input
	///@param[out] a pointor to an array that stores the attributes
	Int getAttr(Int s, Doub * a) {
    takeSamp(attr[s - 1], a, dd);
    return SUCC;
  }

  /**
	 * @brief Get a set of indices corresponding to the set of inputs which
	 * connected with a specific input
	 * @param[in] s Index of a specific input
	 * @param[out] c A list of indices corresponding to the connected inputs
	 * @return Number of the connected inputs 
	 */
  Int getNeighbour(Int s, Vint &c) {
    Int num_c = ne[s - 1][0];
    c.resize(num_c);
    Int cur = 0;
    for (Int i = 0; i < num_c; i++) {
      c[cur] = ne[s - 1][i + 1];
      cur++;
    }
    return cur;
  }
};
/**
 * @class pgpr_data
 * @brief This class provides functionalities of preparing train data,
 * test set and support set
 **/
class pgpr_data
{
private:
  t_state *gst; // Pointer to the status (e.g., selected or not) of domain data
  pgpr_fgp *gp; // Pointer to a GP predictor for active selection of subset
  Int blk_num;	// Number of blocks
  Int samp_num;	// Number of samples 
  FILE *logfp;	// Log file

public:
	/// Pointer to a domain
  pgpr_domain *domain;
	/// Matrix of indices used as training data 
  Mint m_sample;

	/**
	 * @brief Constructor
	 * @param[in] cfg Configure of domain
	 * @param[in] an	Number of blocks
	 * @param[in] tn	Number of samples that are used for training
	 **/
  pgpr_data(pgpr_parse cfg, Int an, Int tn): blk_num(an), samp_num(tn) {
    m_sample.resize(samp_num+1, an);
    domain = new pgpr_domain(cfg.m_connect, cfg.m_xyz);
    gp = new pgpr_fgp(cfg.v_hyp, (Int)cfg.v_param[INP]);
    gst = new t_state(domain->ds);
   
		// As the active selection tends to be long, it's better that the
		// progress are recorded in a log.
    logfp = fopen("pgpr_data.log", "w");
    if (NULL == logfp) {
      throw("failed to open log file\n");
    }
    setlinebuf(logfp);
  
	}

	/**
	 * @brief Get a random subset of data for testing purpose
	 * @param[in] tsize Size of test set
	 * @param[out] test Set of data for test; each row represents a test point.
	 **/
  Int getTestSet(Int tsize, Mdoub &test) {
    test.resize(tsize, domain->dd);
    for (Int i = 0; i < tsize; i++) {
      Int pos  = domain->selRand(gst);
      if(pos < 0) {
        throw("getTestSet: failed to determine test set \n");
      }
      domain->getAttr(pos, test[i]);
    }
    return SUCC;
  }

  /**
	 * @brief Get a set of training data
   * @param[in]  mode different selection algorithms: ALGO1 - Simple
	 * clustering algorithm; ALGO2 - Random walk algorithm
	 **/
  Int getTrainSet(Int mode) {
    switch(mode) {
    case ALGO1:
      getRandomBlk();
      break;
    case ALGO2:
      getRandomWalk();
        break;
    default:
      throw("The method for selecting observations is invalid\n");
    }
    return SUCC;
  }
 
	/**
	 * @brief Get a set of training data based on a simple clustering
	 * scheme 
   * @details This function follows steps (see Definition 5, remark 2):
   * 1) randomly select a set of K central points;
	 * 2) randomly select an unobserved point;
   * 3) allocate it to the closest non-full bin;
   * 4) loop until all bins are full.
	 * @todo the training data is stored in m_sample_{i,j} where i is the
	 * i-th observation and j is the j-th bin. I am considering move
	 * m_sample out of member attribute
	**/
  Int getRandomBlk() {
    pmsg(LEV_DBG, stdout, "Blocking algorithm: simple random clustering\n");
    //size of bins
    Vint bs(blk_num, 0);
    //1.central points
    for (Int i = 0; i < blk_num; i++) {
      Int pos  = domain->selRand(gst);
      m_sample[0][i] = pos;
      bs[i]++;
    }
    //2. observations
    Vdoub dist(blk_num, .0);
    Vint  ai(blk_num);
    for(Int i = 0; i < (samp_num - 1)*blk_num; i++) {
      Int ci  = domain->selRand(gst);
      for (Int k = 0; k < blk_num; k++) {
        dist[k] = domain->getDist(ci, m_sample[0][k]);
        ai[k] = k;
      }
      //Sort in incr order according to distance to centroid.
      bubble_sort(dist, ai);
      //assign point to a non-full bin
      for (Int k = 0; k < blk_num; k++) {
        Int ca = ai[k];
        if (bs[ca] < samp_num) {
          m_sample[bs[ca]][ca] = ci;
          bs[ca]++;
          break;
        }
      }
    }

    return SUCC;
  }
  
	/**
	 * @brief Get a set of training data based on random walk. Note that,
	 * it is applicable only if the domain contains topology information. 
	 * @todo the training data is stored in m_sample_{i,j} where i is the
	 * i-th observation and j is the j-th bin. I am considering move
	 * m_sample out of member attribute
	 **/
  Int getRandomWalk() {
    // 1. randomly locate
    for (Int i = 0; i < blk_num; i++) {
      m_sample[0][i] = domain->selRand(gst);
      gst->upSt(m_sample[0][i]);
    }
    Mdoub cset;
    Vdoub util;
    // 2. loop
    for(Int i = 1; i < samp_num; i++) {
      pmsg(LEV_DBG, logfp, "##########:%d\n", i + 1);
      for (Int j = 0; j < blk_num; j++) {
        Int prev = m_sample[i - 1][j];
        Int cnum = 0;
        //index of candidates
        Vint ci;
        cnum = domain->getNeighbour(prev, ci);
        util.assign(cnum, -1);
        //decide next according to utility vector
        Int mi = RANI(cnum);
        m_sample[i][j] = ci[mi];
      }
    }
    return SUCC;
  }
	
	/**
	 * @brief Actively select an informative subset of points 
	 * @param[in] dset universal set of points 
	 * @param[in] dnum size of the universal set
	 * @param[in] anum size of the subset to be selected 
	 * @param[out] aset the informative subset
	**/
  Int selMaxVar(Mdoub dset, Int dnum, Int anum, Mdoub &aset) {
    Int ddim = dset.ncols();
    Vint flag(dnum, 0);
    Int pt;
    //1. initialized as empty set 
    if (aset.nrows() < anum) {
      aset.resize(anum, ddim);
    }
    pt = RANI(dnum);
    flag[pt] = 1;
    takeSamp(dset[pt], aset[0], ddim);
		//2. actively select a point with largest predictive variance 
		//until the size of set reached a predefined number
    for (Int i = 1; i < anum; i++) {
      //take the complement set
      Mdoub rset(dnum - i, ddim);
      Vint rsi(dnum - i, -1);
      Int ri = 0;
      for (Int j = 0; j < dnum; j++) {
        if(flag[j] == 0) {
          takeSamp(dset[j], rset[ri], ddim);
          rsi[ri] = j;
          ri++;
        }
      }
      //select the largest variance
      Vdoub pmu(dnum - i);
      Vdoub pvar(dnum - i);
      gp->full_reg(aset, i, rset, rset.nrows(), pmu, pvar);
      pt = rsi[argmaxi(pvar)];
      flag[pt] = 1;
      takeSamp(dset[pt], aset[i], ddim);
    }
    return SUCC;
  }

};
#endif

Introduction                         {#mainpage}
========================================================================
The implementation of parallel Gaussian process (GP) regression is based on the following publication:

1. Jie Chen, Nannan Cao, Kian Hsiang Low, Ruofei Ouyang, Colin Keng-Yan Tan & Patrick Jaillet. Parallel Gaussian Process Regression with Low-Rank Covariance Matrix Approximations. In Proceedings of the 29th Conference on Uncertainty in Artificial Intelligence (UAI 2013), Bellevue, WA, Jul 11-15, 2013. It can be found at http://arxiv.org/abs/1305.5826 or http://www.comp.nus.edu.sg/~lowkh/pubs/uai2013.pdf.

2. Kian Hsiang Low, Jiangbo Yu, Jie Chen and Patrick Jaillet. Parallel Gaussian Process Regression for Big Data: Low-Rank Representation Meets Markov Approximation. In Proceedings of the 29th AAAI Conference on Artificial Intelligence (AAAI-15), Austin, TX, Jan 25-29, 2015. It can be found at http://www.comp.nus.edu.sg/~lowkh/pubs/aaai2015.pdf or http://arxiv.org/abs/1411.4510.
========================================================================
1. System requirements 
------------------------------------------------------------------------
Linux/Unix/MacOS X environment

64-bit processors

GNU GCC (4.2.1 or above)

MPICH 3.0.4 (http://www.mpich.org/); we also test in version 1.5

Eigen (http://eigen.tuxfamily.org) for PLMA GP regression, Eigen (version 3.2.3) is already included under lib/

2. Compile  
------------------------------------------------------------------------
To compile all the applications, enter the command:

    make all

All applications will be automatically generated in folder demo.  In
addition, it's also supported to compile each individual application.
For example, to compile the application that prepares data (training
data, test data etc.) for experiments
  
   make prep

For the demonstration of different Gaussian process regression (GPR),
you can use command

    make fgp

to compile the application that demonstrates full Gaussian process
regression;

    make pitc 

to compile the application that demonstrates PITC GP regression; 

    make ppitc 

to compile the application that demonstrates parallel PITC GP
regression; 

    make pic

to compile the application that demonstrates PIC GP regression; 

    make ppic
	
to compile the application that demonstrates PLMA GP regression, 
  make plma

to compile the application that demonstrates PICF-based GP regression.
To clean the compilation environment, use the command:

    make clean



3. Demonstrations
------------------------------------------------------------------------
A bash script can be used to run all applications, using command

    cd demo && bash bat_demo.sh

Basically, the script first prepares all necessary files (training data,
    test data, support set, and hyperparameter file) for experiments;
Then, different GPR algorithms are run sequentially and output the
results (i.e., incurred time, root mean square error (RMSE) and mean
    negative log probability (MNLP) ). For more information about the
arguments of applications, please refer to the comments in the bash
script.

Note when setting the bandwidth to be 0, PLMA would be the same as PPIC.
The results may be different because of different clustering algorithms.
When setting the bandwidth to be the maximum vale (blk - 1), PLMA would be
equal to the FGP.


4. Documentation
------------------------------------------------------------------------
To compile the documentation, the documentation generation tool doxygen 
(http://www.doxygen.org) needs to be installed. Then, enter the home 
directory of our source code and run the command	
	
	make doc

You can refer to the html version documentation by

	cd doc/html 

and use any browser to open index.html; In addition, the pdf version can be 
accessed by

	cd doc/latex && make

and use any pdf viewer to open refman.pdf.



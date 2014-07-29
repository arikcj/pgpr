Introduction                         {#mainpage}
========================================================================
1. System requirements 
------------------------------------------------------------------------
Linux/Unix environment
64-bit processors
GNU GCC (4.2.1 or above)
MPICH 3.0.4 (http://www.mpich.org/); we also test in version 1.5


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



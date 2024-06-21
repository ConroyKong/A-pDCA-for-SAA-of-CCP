# A-pDCA-for-SAA-of-CCP
This code is the implementation of the algorithms and experiments in the paper "A Proximal DC Algorithm for Sample Average Approximation of Chance Constrained Programming".

* `Portfolio`: the source code for the experiments of "VaR-Constrained Portfolio Selection Problem".
* `PTP`: the source code for the experiments of "Probabilistic Transportation Problem with Convex Objective".
* `CCSCP`: the source code for the experiments of "Probabilistic Transportation Problem with Non-Convex Objective".
* `Norm Optimization Problem`: the source code for the experiments of "Linear Optimization with Nonlinear Chance Constraint".

### Introduction
This code can be run in Matlab R2022b and Gurobi 9.5.2. 

* `main_xxxx.m`: runing file.
* `MIP.m`: the implementation of the mixed-integer program (MIP) in paper "Solving chance-constrained stochastic programs via sampling and integer programming".
* `CVaR.m`: the implementation of the CVaR in paper "Convex approximations of chance constrained programs".
* `BiCVaR.m`: the implementation of the bisection-based CVaR method in paper "An augmented Lagrangian decomposition method for chance-constrained optimization problems".
* `DCA.m`: the implementation of Algorithm 1 without proximal term in our paper.
* `pDCA.m`: the implementation of Algorithm 1 in our paper.
* `SCA.m`:  the implementation of the successive convex approximation method (SCA) in paper "Squential convex approximations to joint chance constrained programs: A Monte Carlo approach".
* `ALDM.m`, `ALDM_update_x.m`, `ALDM_update_y.m`: the implementation of the augmented Lagrangian decomposition method (ALDM) in paper "An augmented Lagrangian decomposition method for chance-constrained optimization problems".
* `post_processing.m`: the code for calculating the maximum, minimum and mean values of the indicators.
* `risk_level.m`: the code for calculating the risk level of each method.
* `Norm Optimization Problem/gensample.m`: the code for generating the $d \times m$ matrix of random variables $\xi$.

### How to get the results
* To run the experiments of "VaR-Constrained Portfolio Selection Problem", please `main_portfolio.m`.
* To run the experiments of "Probabilistic Transportation Problem with Convex Objective", please `main_PTP.m`
* To run the experiments of "Probabilistic Transportation Problem with Non-Convex Objective", please `main_CCSCP.m`.
* To run the experiments of "Linear Optimization with Nonlinear Chance Constraint", please `main_NormOpt.m`.


### Citation
If you want to reference our paper in your research, please consider citing us by using the following BibTeX:

```bib
@Article{Wang2023,
  author        = {Peng Wang and Rujun Jiang and Qingyuan Kong and Laura Balzano},
  title         = {Difference-of-Convex Reformulation for Chance Constrained Programs},
  year          = {2023},
  month         = jan,
  archiveprefix = {arXiv},
  eprint        = {2301.00423},
  primaryclass  = {math.OC},
}
```

If you have any questions, please contact us.

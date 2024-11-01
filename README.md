# master-thesis

pdf and codes for my master thesis

CODE FOLDER

main codes:
- multiple: run the simulations, merge the results, create and save the tables
- roots_plot: plot the VAR and VMA roots in different scenarios

functions:
- sim: run simulations for different scenarios
- coeff_loop: choose the right cell of VAR parameters
- A2C: convert the cell of VAR parameters in the companion matrix
- VARMA_IRF: compute the IRF of a VARMA process
- loop: extract the C.I. and the IRF estimates using LP and lag-augmented LP for all persistence levels, given all the other parameters
- LP_loop: simulate a DGP given the parameters, compute LP (classical or lag-augmented) and the C.I. using asymptotics, recursive and wild bootstrap
- VARMA_sim: simulate a VARMA process
- LP_est: estimate the LP and its s.e. using OLS and Newey-West s.e. (classical heteroskedasticity robust s.e. are employed when lag-augmented LP is used)
- VAR_est: estimate a VAR and its s.e. by OLS (it does not use HC s.e.)
- VAR_boot: regenerate recursively B bootstrapped samples (either by iid or wild bootstrap). The bootstrapped series are returned in a single array in which the first M columns are the first bootstrap sample and so on
- VARGARCH_sim: simulate a VAR with MGARCH errors
- GARCH_sim: simulate a MGARCH process
- VARH_sim: simulate a VAR with structural break in the errors covariance matrix
- tab: converts the C.I. arrays in coverage tables and creates csv and tex files
- cov2tab: converts the logical arrays ("true value is in the bounds?") into coverage tables
- VAR_roots: compute the VAR roots given the VAR params cell
- VMA_roots: compute the VMA roots given the VMA params cell

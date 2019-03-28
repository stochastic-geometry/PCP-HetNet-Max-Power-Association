# Matlab code for the computation of coverage probability in  the unified HetNet model with Poisson point process and Poisson cluster process
---
This repository contains the matlab scripts used to generate the results of the paper [https://arxiv.org/pdf/1812.01830.pdf](https://arxiv.org/pdf/1812.01830.pdf) presented in Section-IV. 

For the numerical evaluation, we choose a two tier network where BSs in tier-1 are distributed as Poisson cluster process (PCP) and BSs in tier 2 are distributed as Poisson point process (PPP). The coverage probability (P_c) is evaluated for different values of coverage threshold from simulation and analysis. 

Run 'CoverageSimulation.m' to generate coverage probability by Monte Carlo simulation of the network. 

Run ComputeCoverage[PCP]Type[X]user.m to generate coverage probability by evaluating the theoritical expressions derived in the paper (see Theorems 1 and 2), where PCP = TCP or MCP, X = 1 or 2.

Email to csaha@vt.edu for further questions/issues.

Please use the following bibtex code to cite the paper if parts of the scripts are reused. 

```
@article{saha2018unified,
  title={Unified Analysis of HetNets using Poisson Cluster Process under Max-Power Association},
  author={Saha, Chiranjib and Dhillon, Harpreet S and Miyoshi, Naoto and Andrews, Jeffrey G},
  note={available online: arXiv/abs/1812.01830},
  year={2018}
}
```




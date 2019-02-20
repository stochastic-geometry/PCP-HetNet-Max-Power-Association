function [pc1,pc2] = ComputeCoverageThomasType1user(P,tau_dB,Lambda,m,sigma)
% =======================================================================================
% Author: Chiranjib Saha, Harpreet S. Dhillon
% The code is a part of the repository 'PCP-HetNet-Max-Power-Association'
% which is the collection of all Matlab codes used to generate the results
% of the following paper:
%@article{saha2018unified,
%  title={Unified Analysis of HetNets using Poisson Cluster Process under Max-Power Association},
%  author={Saha, Chiranjib and Dhillon, Harpreet S and Miyoshi, Naoto and Andrews, Jeffrey G},
%  note={available online: arXiv/abs/1812.01830},
%  year={2018}
%}
% =======================================================================================
% This function computes coverage probability of two tier HetNet described in  for Type-1
% users when Tier 1 is modelled as a Thomas cluster process.
% Refer to Section-IV for the details of the system model. 
% Input ---
% P: 1x2 vector specifying the power levels
% tau_dB: Coverage theshold in dB
% Lambda: 1x2 vector specifying BS intensities
% m: average number of BSs of the TCP
% sigma: cluster standard deviation of TCP
% Output: per tier coverage
P_1  = P(1);
P_2 =  P(2);


l_p_2 = Lambda(2);
l_p_1 = Lambda(1);

alpha = 4; 
tau = 10^(tau_dB/10);

f_d1 = @(x,r)    ricepdf(x,r,sigma);
%% When user connects to tier 1 %% 
 %%% P-s %%
 P_11 = 1;
 P_21 = (P_2/P_1)^(1/alpha);
 %%% C-s %%%(FOR SUM-Product)
 C1 = @(r,s) exp(-m*(1-(integral(@(u)f_d1(u,s)./(1+tau*(P_11*r./u).^alpha),P_11*r,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));
 C_arr1 = @(r,s) arrayfun(@(r,s) C1(r,s),r,s);
 %%% T-s %%%%
 T1  = @(r) 2*pi*l_p_1.*integral(@(s)f_d1(r,s).*C_arr1(r,s).*s,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3);
 %%% M-s: combine all PGFLs %%%%%
 M1 = @(r)  exp(-2*pi*l_p_1*integral(@(s) ...
                          (1-C_arr1(r,s)).*s,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
 M2 = @(r)exp(-pi*l_p_2*P_21^2.*r.^2.*(1+2*tau./(alpha-2).*hypergeom([1,1-2/alpha],2-2/alpha,-tau)));
 f1 = @(r) T1(r).*M1(r).*M2(r);
 f1_arr = @(r) arrayfun(@(r)f1(r),r);
 pc1 = m*integral(@(r) f1_arr(r),0,20*sigma,'reltol',1e-3,'abstol',1e-3);
%% When user connects to tier 2 %% 
 %%% P-s %%
 P_12 = (P_1/P_2)^(1/alpha);
 P_22 = 1;
 %%% C-s %%%(FOR SUM-Product)
 C1 = @(r,s) exp(-m*(1-(integral(@(u)f_d1(u,s)./(1+tau*(P_12*r./u).^alpha),P_12*r,40*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))));
 C_arr1 = @(r,s) arrayfun(@(r,s) C1(r,s),r,s);
 %%% M-s: combine all PGFLs %%%%%
 M1 = @(r)  exp(-2*pi*l_p_1*integral(@(s) ...
                          (1-C_arr1(r,s)).*s,0,20*sigma,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));
 M2 = @(r)2*pi*l_p_2.*r.*exp(-pi*l_p_2*P_22^2.*r.^2.*(1+2*tau./(alpha-2).*hypergeom([1,1-2/alpha],2-2/alpha,-tau)));
 f2 = @(r) M1(r).*M2(r);
 f2_arr = @(r) arrayfun(@(r)f2(r),r);
 pc2 = integral(@(r) f2_arr(r),0,20*sigma,'reltol',1e-3,'abstol',1e-3);
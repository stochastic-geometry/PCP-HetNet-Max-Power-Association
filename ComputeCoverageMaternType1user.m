function [pc1,pc2] = ComputeCoverageMaternType1user(P,tau_dB,Lambda,m,R)
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
% users when Tier 1 is modelled as a Matern cluster process.
% Refer to Section-IV for the details of the system model. 
% Input ---
% P: 1x2 vector specifying the power levels
% tau_dB: Coverage theshold in dB
% Lambda: 1x2 vector specifying BS intensities
% m: average number of BSs of the TCP
% R: cluster radius of MCP
% Output: per tier coverage. Coverage probability can be computed as Pc= pc1+pc2.
P_1  = P(1);
P_2 =  P(2);


l_p_2 = Lambda(2);
l_p_1 = Lambda(1);

alpha = 4; 
tau = 10^(tau_dB/10);

%%%%%
  % if  (x<r-z)&&(z<=r)
 chi_1 = @(x)  2.*x/R^2;
% if (x>abs(r-z))&&(x<r+z)
 chi_2 = @(x,z)  2.*x/(pi*R^2).*acos((x.^2+z.^2-R^2)./(2.*x.*z));
 chi_3= chi_2;
 %% 1st Part %%
 % if z<r, x<r-z %
 C_1 = @(x,z) exp(- m*(1-integral(@(u)chi_1(u)./(1+tau*(x./u)^alpha),x,R-z,'arrayvalued',true)...
    -  integral(@(u)chi_2(u,z)./(1+tau*(x./u)^alpha),R-z,R+z,'arrayvalued',true)));

 % if z<r, r-z<x<z+r %
 C_2 = @(x,z)  exp(-m*(1-integral(@(u)chi_2(u,z)./(1+tau*(x./u)^alpha),x,R+z,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)));
 % if z>r, 0<x<z+r %
 C_3 = @(x,z)  exp(-m*(1-integral(@(u)chi_2(u,z)./(1+tau*(x./u)^alpha),max(x,z-R),R+z,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))); 
 % if x>r+z , irrespective of z
 C_4 = @(x,z) exp(-m);
 %% for sanity check of C 
 % indicates which piece of C will be active given an ordered pair (x,z)
 indicator_C = @(x,z) 1.*(z<R).*(x<R-z)+ 2.*(z<R).*(x>R-z).*(x<z+R)+   3.*(z>R).*(x<z+R)   +  4*(x>R+z);
 %% Vectorizing C 
 C_1arr= @(x,z)arrayfun(@(x,z)C_1(x,z),x,z);
 C_2arr= @(x,z)arrayfun(@(x,z)C_2(x,z),x,z);
 C_3arr= @(x,z)arrayfun(@(x,z)C_3(x,z),x,z);
 C_4arr=  C_4;
 %% the T-function
 % x<r 
 T_1 = @(x) 2*pi*l_p_1 * (  chi_1(x).*integral(@(z) z.* C_1arr(x,z),0,R-x,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) z.*chi_2(x,z).* C_2arr(x,z),R-x,R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                           integral(@(z) z.*chi_3(x,z).* C_3arr(x,z),R,R+x,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)); 
 % r<x<2r
 T_2 = @(x) 2*pi*l_p_1 * ( integral(@(z) z.*chi_2(x,z).* C_2arr(x,z),x-R,R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
        integral(@(z) z.*chi_3(x,z).* C_3arr(x,z),R,x+R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)); 
% x>2r
 T_3 = @(x) 2*pi*l_p_1 * ( integral(@(z) z.*chi_3(x,z).* C_3arr(x,z),x-R,x+R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3));%+...
 
 %% the M-function
 % x<r
 M_1 = @(x)  exp(-2*pi*l_p_1*...
             ( integral(@(z) ...
                          (1-C_1arr(x,z)).*z,0,R-x,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       +integral(@(z) ...
                          (1-C_2arr(x,z)).*z,R-x,R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                      + integral(@(z)(1-C_3arr(x,z)).*z,R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...   
                     )  );
% r<x<2*r
 M_2 = @(x)  exp(-2*pi*l_p_1*...
              (integral(@(z) ...
                          (1-C_4arr(x,z)).*z,0,x-R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_2arr(x,z)).*z,x-R,R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_3arr(x,z)).*z,R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       ));
% x>2*r 
 M_3 = @(x)  exp(-2*pi*l_p_1*...
              (integral(@(z) ...
                          (1-C_4arr(x,z)).*z,0,x-R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_3arr(x,z)).*z,x-R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       ));
 M_PPP = @(r)exp(-pi*l_p_2*P_21^2.*r.^2.*(1+2*tau./(alpha-2).*hypergeom([1,1-2/alpha],2-2/alpha,-tau)));
%%
fnc1 = @(r) T_1(r).*M_1(r).*M_PPP(r);
fnc2 = @(r) T_2(r).*M_2(r).*M_PPP(r);
fnc3 = @(r) T_3(r).*M_3(r).*M_PPP(r);
fnc1_arr = @(r) arrayfun(@(r)fnc1(r),r);
fnc2_arr = @(r) arrayfun(@(r)fnc2(r),r);
fnc3_arr = @(r) arrayfun(@(r)fnc3(r),r);

pc_11 = m*integral(@(r)fnc1_arr(r),0,R,'reltol',1e-3,'abstol',1e-3);
pc_12 = m*integral(@(r)fnc2_arr(r),R,2*R,'reltol',1e-3,'abstol',1e-3);
pc_13 = m*integral(@(r)fnc3_arr(r),2*R,20*R,'reltol',1e-3,'abstol',1e-3);
pc1 = pc_11+pc_12+pc_13;
%%%%%
P_12 = (P_1/P_2)^(1/alpha);
P_22 = 1;
  % if z<R, x<P_12^-1*(R-z) %
C_1 = @(x,z) exp(- m*(1-integral(@(u)chi_1(u)./(1+tau*(P_12*x./u)^alpha),P_12*x,R-z,'arrayvalued',true)...
    -  integral(@(u)chi_2(u,z)./(1+tau*(P_12*x./u)^alpha),R-z,R+z,'arrayvalued',true)));
  % if z<R, P_12^-1*(R-z)<x<P_12^-1*(z+r) %
C_2 = @(x,z)  exp(-m*(1-integral(@(u)chi_2(u,z)./(1+tau*(P_12*x./u)^alpha),P_12*x,R+z,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)));
  % if z>R, 0<x<P_12^-1*(z+R) %
C_3 = @(x,z)  exp(-m*(1-integral(@(u)chi_2(u,z)./(1+tau*(P_12*x./u)^alpha),max(P_12*x,z-R),R+z,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3))); 
  % if x>P_12^-1*(R+z) , irrespective of z
C_4 = @(x,z) exp(-m);
C_1arr= @(x,z)arrayfun(@(x,z)C_1(x,z),x,z);
C_2arr= @(x,z)arrayfun(@(x,z)C_2(x,z),x,z);
C_3arr= @(x,z)arrayfun(@(x,z)C_3(x,z),x,z);
C_4arr=  C_4;
  %x<P_12^-1*R
M_1 = @(x)  exp(-2*pi*l_p_1*...
             ( integral(@(z) ...
                          (1-C_1arr(x,z)).*z,0,(R-P_12*x),'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       +integral(@(z) ...
                          (1-C_2arr(x,z)).*z,(R-P_12*x),R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                      + integral(@(z)...
                          (1-C_3arr(x,z)).*z,R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...   
                     )  );
  %P_12^-1*R<x<2*P_12^-1*R
M_2 = @(x)  exp(-2*pi*l_p_1*...
              (integral(@(z) ...
                          (1-C_4arr(x,z)).*z,0,P_12*x-R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_2arr(x,z)).*z,P_12*x-R,R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_3arr(x,z)).*z,R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       ));
  %x>2*P_12^-1*R
M_3 = @(x)  exp(-2*pi*l_p_1*...
              (integral(@(z) ...
                          (1-C_4arr(x,z)).*z,0,P_12*x-R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)+...
                          integral(@(z) ...
                          (1-C_3arr(x,z)).*z,P_12*x-R,20*R,'arrayvalued',true,'reltol',1e-3,'abstol',1e-3)...
                       ));
M_PPP = @(r)2*pi*l_p_2.*r.*exp(-pi*l_p_2*P_22^2.*r.^2.*(1+2*tau./(alpha-2).*hypergeom([1,1-2/alpha],2-2/alpha,-tau)));
  fnc1 = @(r) M_1(r).*M_PPP(r);
  fnc2 = @(r) M_2(r).*M_PPP(r);
  fnc3 = @(r) M_3(r).*M_PPP(r);
  fnc1_arr = @(r) arrayfun(@(r)fnc1(r),r);
  fnc2_arr = @(r) arrayfun(@(r)fnc2(r),r);
  fnc3_arr = @(r) arrayfun(@(r)fnc3(r),r);
   pc_11 = integral(@(r)fnc1_arr(r),0,P_12^-1*R,'reltol',1e-3,'abstol',1e-3);
   pc_12 = integral(@(r)fnc2_arr(r),P_12^-1*R,2*P_12^-1*R,'reltol',1e-3,'abstol',1e-3);
   pc_13 = integral(@(r)fnc3_arr(r),2*P_12^-1*R,20*P_12^-1*R,'reltol',1e-3,'abstol',1e-3);

pc2 = pc_11+pc_12+pc_13;
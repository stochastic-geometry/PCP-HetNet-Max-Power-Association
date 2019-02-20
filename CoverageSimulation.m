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
% This code simulates the  Fig. 2-a (indicated as Circles) 
% Set PCP_type equal to 'TCP' or 'MCP'
% Set user_type = 1 or 2 according to type-I or type-II users.
clear all;close all;


lambda_SBS_t1= 2.5e-5; %/km^2 
lambda_SBS_t2= 1e-6; %/km^2 

m_1=4;  % average number of points per cluster for the MCP
%PCP_type = 'TCP';
%sigma = 20;% cluster std

PCP_type = 'MCP';
R = 20;% cluster radius

user_type = 1;

% trasnmit powers
P_1  = 1;
P_2 = 1000;

betaConst=4; %path-loss exponent  
%Simulation section
simNumb=1e5;%1e6; %number of simulations
diskRadius=3500; %km^2 %radius of simulation disk region (has to be larger when fading is incorporated)
%diskRadius=10^3; %radius of simulation disk region (has to be larger when fading is incorporated)
tier_count=0;
dist_count=0;
served_by_other_point_covered_count=0;
diskArea=pi*diskRadius^2;

%power_ratio= 10^-4.0; 
%W = 10e6;
tau_all =  -20:4:20;
count_tau = 0;
for count_tau=1:length(tau_all)
 rho = tau_all(count_tau);%50e3;
 cov_sbs_sim_t1 = 0 ; 
 cov_sbs_sim_t2 = 0 ;
 for count_sim=1:simNumb
  %  count_sim
    randNumb_SBS_t1=poissrnd(lambda_SBS_t1*diskArea);
    randNumb_SBS_t2=poissrnd(lambda_SBS_t2*diskArea);
 %% TIER 1
  %% Generate BS
   theta = rand(randNumb_SBS_t1,1)*(2*pi);
   r = diskRadius*sqrt(rand(randNumb_SBS_t1,1));
   x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y_1 =  r.*sin(theta);   %%%************************************************
   
  if user_type ==2   
   if strcmp(PCP_type,'TCP')
        x_1 =  [x_1; sigma*randn];   %%%*****shifting origin to receiver location*******
        y_1 =  [y_1; sigma*randn];
   elseif strcmp(PCP_type,'MCP')
       theta1 = rand*(2*pi);
       r1 = R*sqrt(rand);
       x_1 = [x_1;r1*cos(theta1)];
       y_1 = [y_1;r1*sin(theta1)];
   else
     error('Invalid type of PCP. It should be either TCP or MCP');
   end
   randNumb_SBS_t1 = randNumb_SBS_t1+1;
  end
  SBS_location = [x_1,y_1];
   
   %% this is the typical cluster center %% 
   %x_0= sigma*randn(1,2);
   UE_cc = SBS_location;
  % SBS_location = UE_cc;
   UE_location_all=[];
   total_user_count=0;
   no_users_t1= poissrnd(m_1,randNumb_SBS_t1,1);
   %no_users(1) = no_users(1)+1; %% adding the count of the typical user 
   r = no_users_t1;
   x = UE_cc;
   t = r > 0;
   a = cumsum(r(t));
   b = zeros(1,a(end));
   b(a - r(t) + 1) = 1;
   x1 = UE_cc(t,:);
   c = x1(cumsum(b),:);
   if strcmp(PCP_type,'TCP')
     user_pos=sigma*randn(sum(no_users_t1),2);%R*randn(sum(no_users_t1),2);
   elseif strcmp(PCP_type,'MCP')
     theta1 = rand(sum(no_users_t1),1)*(2*pi);
     r1 = R*sqrt(rand(sum(no_users_t1),1));
     user_pos = [ r1.*cos(theta1) r1.*sin(theta1)];
   else
     error('Invalid type of PCP. It should be either TCP or MCP');
   end
   cc_location_rep = x1(cumsum(b),:);
   BS_location_all_t1=cc_location_rep+ user_pos;
 
  %% TIER 2
  %% Generate BS
   theta = rand(randNumb_SBS_t2,1)*(2*pi);
   r = diskRadius*sqrt(rand(randNumb_SBS_t2,1));
   x_1 =  r.*cos(theta);   %%%*****shifting origin to receiver location*******
   y_1 =  r.*sin(theta);   %%%************************************************
   SBS_location = [x_1,y_1];
   no_users_t2 = randNumb_SBS_t2;
   BS_location_all_t2= SBS_location;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   macro_dist_t1 = sqrt(BS_location_all_t1(:,1).^2+BS_location_all_t1(:,2).^2);
   macro_dist_t2 = sqrt(BS_location_all_t2(:,1).^2+BS_location_all_t2(:,2).^2);
   h_m_t1 = exprnd(1,sum(no_users_t1),1);
   h_m_t2 = exprnd(1,sum(no_users_t2),1);
   [mindist_t1, minindex_t1] = min(macro_dist_t1);
   [mindist_t2, minindex_t2] = min(macro_dist_t2);
   I_t1 = sum(P_1*h_m_t1.* macro_dist_t1.^-betaConst);
   I_t2 = sum(P_2*h_m_t2.* macro_dist_t2.^-betaConst);
   if P_1*mindist_t1^-betaConst>P_2*mindist_t2^-betaConst
       association_event = 1;
       S_bs = P_1* h_m_t1(minindex_t1).*mindist_t1^-betaConst;
       %I_t2 = 0 ; 
       SIR =  S_bs / (I_t1+I_t2-S_bs);
       if SIR>10^(rho/10)
         cov_sbs_sim_t1 = cov_sbs_sim_t1+1;
      end
     %cov_sbs_sim_t1/count_sim
   else 
       association_event = 2;
       S_bs = P_2* h_m_t2(minindex_t2).*mindist_t2^-betaConst;
       %I_t2 = S_bs ;
       SIR =  S_bs / (I_t1+I_t2-S_bs);
       if SIR>10^(rho/10)
         cov_sbs_sim_t2 = cov_sbs_sim_t2+1;
       end
     % cov_sbs_sim_t2/count_sim
   end
 %  fprintf('\n Coverage: tier 1 = %f, tier 2 = %f\n',cov_sbs_sim_t1/count_sim,cov_sbs_sim_t2/count_sim);
end 
   Pc_sim(count_tau) =( cov_sbs_sim_t1+ cov_sbs_sim_t2)/count_sim;
end
% save FileName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo of the numerical experiment in the paper "Joint Community Detection 
% and Rotational Synchronization via Semidefinite Programming." by Yifeng
% Fan, Yuehaw Khoo, and Zhizhen Zhao.
%
% Currently, we only support SO(2) and SO(3) transformations, i.e. d = 2 or
% d = 3.
% 
% We use MOSEK as the SDP solver (https://www.mosek.com/), which is
% included in the CVX package (http://cvxr.com/cvx/).
%
% Yifeng Fan (yifengf2@illinois.edu), Dec 2021

clear 
addpath(genpath('./'))
rng('default')

%%% Parameter setting
m_list = [50, 50]; % The cluster sizes
d = 2; 
p = 1;
q = 0;
K = numel(m_list); 

%%% Generate the observation matrix A
[A, ~, M] = gen_observation(m_list, p, q, d); 

%%% Solve by SDP
% M_SDP = sync_SDP_equal(A, d); % Use this when two cluster sizes are equal
% M_SDP = sync_SDP_unequal(A, d, m_list); % Use this when two cluster sizes are unequal
M_SDP = sync_SDP_unknown(A, d); % Use this when two cluster sizes are unknown

%%% Check the recovery error 
% error_SDP_known(M, M_SDP) % Use this when two cluster sizes are known
error_SDP_unknown(M, M_SDP, m_list, d) % Use this when two cluster sizes are unknown
    
    


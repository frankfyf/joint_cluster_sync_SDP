%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SDP for two-unequal sized clusters
%
% Inputs:
%      A: The observation matrix, should be of size nd by nd
%      d: the dimension of rotational transformation in SO(d)
%      m_list: The list of cluster sizes, should be of length 2  
% Outputs:
%      M: The SDP solution

function [ M_SDP ] = sync_SDP_unequal( A, d, m_list )

%Parameters
n = size(A,1)/d;
m = max(m_list); % The maximum cluster size

% SDP solving by MOSEK in CVX
cvx_solver  
cvx_begin sdp 
    variable M_SDP(n*d,n*d) semidefinite;
    expression M_fro(n,n)
    for i = 1:n
        for j = 1:n
            M_fro(i,j) = norm(M_SDP((i-1)*d+1:i*d, (j-1)*d+1:j*d), 'fro');
        end
    end
    maximize( trace(A*M_SDP) );
    subject to
         for i = 1:n
            M_SDP((i-1)*d+1:i*d, (i-1)*d+1:i*d) == eye(d);
            sum(M_fro(i,:)) <= m*sqrt(d);
         end
         sum(M_fro(:)) <= sum(m_list.^2)*sqrt(d);
cvx_end

end

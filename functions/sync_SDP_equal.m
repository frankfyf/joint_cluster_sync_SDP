%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SDP for two-equal sized clusters
%
% Inputs:
%      A: The observation matrix, should be of size nd by nd
%      d: the dimension of rotational transformation in SO(d)
% Outputs:
%      M: The SDP solution

function [ M_SDP ] = sync_SDP_equal( A, d )

%Parameters
n = size(A,1)/d;

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
            sum(M_fro(i,:)) <= n*sqrt(d)/2;
         end
cvx_end

end


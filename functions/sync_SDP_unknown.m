%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SDP for two clusters with unknown cluster sizes
%
% Inputs:
%      A: The observation matrix, should be of size nd by nd
%      d: the dimension of rotational transformation in SO(d)
% Outputs:
%      M: The SDP solution

function [ M_SDP ] = sync_SDP_unknown( A, d )

%Parameters
n = size(A,1)/d;

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
         for i = 1:d
             for j = 1:d
                 if i == j
                     sum(M_SDP(((j-1)*n*d+i):(n*(d^2) + d):end)) == 2;
                 else
                     sum(M_SDP(((j-1)*n*d+i):(n*(d^2) + d):end)) == 0;
                 end
             end
         end
         for i = 1:n
            sum(M_fro(i,:)) <= sqrt(d);
         end
cvx_end

end


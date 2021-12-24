%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the observation matrix A, currently we support SO(2) and SO(3)
% transformations. We assume the first m_1 nodes form the first cluster,
% the following m_2 nodes form the second one, and so on...
%
% Inputs:
%      m_list: The list of cluster sizes, should be of length K, where K is
%              the number of clusters
%      p: The probability of connection within clusters
%      q: The probability of connection across different clusters
%      d: the dimension of rotational transformation in SO(d)
% Outputs:
%      A: The observation matrix, of size nd by nd
%      V: The concatenation of the true rotation of each node, is of size nd by d
%      A_true: The ground truth of A

function [ A, V, A_true ] = gen_observation( m_list, p, q, d )

% Parameters
n = sum(m_list); % The total number of data points
K = numel(m_list); % The number of blocks

% Generate the SO(d) transformations
V = zeros(d*n,d);
if d == 2
    for i = 1:n
        tmp_angle = rand(1)*2*pi;
        V((i-1)*d+1:i*d, :) = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)];
    end
elseif d == 3
    for i = 1:n
        V((i-1)*d+1:i*d,:) = q_to_rot(qrand(1));
    end
else
    error("The dimension d should be 2 or 3")
end

% Generate the observation matrix 
A = cell(K,K);
A_true = zeros(n*d, n*d);
count = 0;
for i = 1:K
    for j = i:K
        if j == i % The block within each community
            
            tmp = V(count+1:count+d*m_list(i),:)*V(count+1:count+d*m_list(i),:)';
            A_true(count+1:count+d*m_list(i), count+1:count+d*m_list(i)) = tmp;
            for l1 = 1:m_list(i)
                for l2 = l1+1:m_list(i)
                    if rand(1) > p
                        tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = zeros(d,d);
                    end
                end
            end
            tmp = triu(tmp);
            tmp = tmp + tmp';
            tmp(1:m_list(i)*d+1:end) = 1;
            A{i,j} = tmp;
            
        else % The block across different communities
            
            tmp = zeros(d*m_list(i), d*m_list(j));
            for l1 = 1:m_list(i)
                for l2 = 1:m_list(j)
                    if rand(1) <= q
                        if d == 2
                            tmp_angle = rand(1)*2*pi;
                            tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = [cos(tmp_angle), sin(tmp_angle); -sin(tmp_angle), cos(tmp_angle)];
                        elseif d == 3
                            tmp_rot = q_to_rot(qrand(1));
                            [ u, ~, v] = svd(tmp_rot);
                            tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = u*v';
                        end
                    else
                        tmp((l1-1)*d+1:l1*d, (l2-1)*d+1:l2*d) = zeros(d,d);
                    end
                end
            end
            A{i,j} = tmp;   
            A{j,i} = tmp';
            
        end
    end
    count = count + d*m_list(i);
end
A = cell2mat(A);
A = sparse(A);
    
end


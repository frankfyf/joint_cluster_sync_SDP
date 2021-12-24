function [ e ] = error_SDP_unknown( M, M_SDP, m_list, d )

    K = numel(m_list);
    count = 0;

    for k = 1:K
        M(d*count+1:d*(count+m_list(k)), d*count+1:d*(count+m_list(k))) = M(d*count+1:d*(count+m_list(k)), d*count+1:d*(count+m_list(k)))/m_list(k);
        count = count + m_list(k);
    end
    
    e = norm(M - M_SDP, 'fro');
end


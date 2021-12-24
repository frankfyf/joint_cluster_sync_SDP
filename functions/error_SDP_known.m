function [ e ] = error_SDP_known( M, M_SDP )
    e = norm(M - M_SDP, 'fro');
end


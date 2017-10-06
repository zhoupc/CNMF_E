function W = compress_weights(weights)
% compress weights into a  sparse matrix 

%% Author: Pengcheng Zhou, Columbia University, 2017 
%% Email: zhoupc1988@gmail.com 

[d1, d2] = size(weights); 
d = d1*d2; 
nmax = 300; 
ii = zeros(d*nmax, 1); 
jj = ii; 
ss = ii; 
k = 0; 
for m=1:d
    temp = weights{m}; 
    tmp_k = size(temp,2); 
    ii(k+(1:tmp_k)) = m; 
    jj(k+(1:tmp_k)) = temp(1,:); 
    ss(k+(1:tmp_k)) = temp(2,:); 
    k = k + tmp_k; 
end 
ii((k+1):end) = []; 
jj((k+1):end) = []; 
ss((k+1):end) = []; 
W = sparse(ii, jj, ss, d, d); 
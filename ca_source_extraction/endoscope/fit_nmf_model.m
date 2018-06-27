function [b, f] = fit_nmf_model(Y, nb, A, C, b_old, f_old, thresh_outlier, sn, ind_patch)
% fit a patched data with NMF
if ~exist('A', 'var') || isempty(A)
    A = ones(d,1); 
    C = zeros(1, T); 
elseif issparse(A)
    A = full(A); 
end
Y = double(Y); 
B = Y - A*C; 
B_old = b_old*f_old; 
tmp_B = B(ind_patch, :); 
ind_outlier = bsxfun(@gt, tmp_B, bsxfun(@plus, B_old, thresh_outlier*reshape(sn, [], 1))); 
tmp_B(ind_outlier) = B_old(ind_outlier);  
B(ind_patch, :) = tmp_B; 
clear tmp_B B_old ind_outlier;

% if (norm(sum(b_old(:)))==0) || (norm(sum(f_old(:)))==0)
[b, f] = nnmf(B, nb);
% else
%     [b, f] = nnmf(B, nb, 'w0', b_old, 'h0', f_old);
% end
b = b(ind_patch, :);

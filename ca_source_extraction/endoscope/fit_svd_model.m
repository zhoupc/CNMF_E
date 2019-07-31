function [b, f, b0] = fit_svd_model(Y, nb, A, C, b_old, f_old, thresh_outlier, sn, ind_patch)
% fit a patched data with SVD
[d, T] = size(Y); 

Ymean = mean(Y,2); 
if ~exist('A', 'var') || isempty(A)
    A = ones(d,1); 
    C = zeros(1, T); 
elseif issparse(A)
    A = full(A); 
end
Cmean = mean(C, 2); 
Y = bsxfun(@minus, double(Y), Ymean); 
C = bsxfun(@minus, C, Cmean);

if ~exist('ind_patch', 'var')
    ind_patch = true(size(A,1), 1); 
end

if nb==0
    b = 0;
    f = 0;
    b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean-b*mean(f, 2);
    return;
end

Bf = Y - A*C;

if exist('thresho_outlier', 'var')
    Bf_old = b_old*f_old;
    tmp_Bf = Bf(ind_patch, :);
    ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1)));
    tmp_Bf(ind_outlier) = Bf_old(ind_outlier);
    Bf(ind_patch, :) = tmp_Bf;
    clear tmp_Bf Bf_old ind_outlier;
end

[u, s, v] = svdsecon(bsxfun(@minus, Bf, mean(Bf, 2)), nb); 
b = u(ind_patch, :)*s; 
f = v'; 
b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean-b*mean(f, 2);


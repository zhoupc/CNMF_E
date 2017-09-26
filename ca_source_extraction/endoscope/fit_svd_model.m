function [b, f, b0] = fit_svd_model(Y, nb, A, C, b_old, f_old, thresh_outlier, sn, ind_patch)
% fit a patched data with SVD
Ymean = mean(Y,2); 
Cmean = mean(C, 2); 
Y = bsxfun(@minus, double(Y), Ymean); 
C = bsxfun(@minus, C, Cmean);
Bf = Y - A*C; 
Bf_old = b_old*f_old; 
tmp_Bf = Bf(ind_patch, :); 
ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1))); 
tmp_Bf(ind_outlier) = Bf_old(ind_outlier);  
Bf(ind_patch, :) = tmp_Bf; 
clear tmp_Bf Bf_old ind_outlier; 

[u, s, v] = svdsecon(bsxfun(@minus, Bf, mean(Bf, 2)), nb); 
b = u(ind_patch, :)*s; 
f = v'; 
b0 = Ymean(ind_patch) - A(ind_patch,:)*Cmean-b*mean(f, 2);


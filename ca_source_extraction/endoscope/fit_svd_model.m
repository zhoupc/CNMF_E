function [b, f, b0] = fit_svd_model(Ypatch, nb, A_patch, C_patch, b_old, f_old, thresh_outlier, sn, ind_patch)
% fit a patched data with SVD
Ymean = mean(Ypatch,2); 
Cmean = mean(C_patch, 2); 
Ypatch = bsxfun(@minus, double(Ypatch), Ymean); 
C_patch = bsxfun(@minus, C_patch, Cmean);
Bf = Ypatch - A_patch*C_patch; 
Bf_old = b_old*f_old; 
tmp_Bf = Bf(ind_patch, :); 
ind_outlier = bsxfun(@gt, tmp_Bf, bsxfun(@plus, Bf_old, thresh_outlier*reshape(sn, [], 1))); 
tmp_Bf(ind_outlier) = Bf_old(ind_outlier);  
Bf(ind_patch, :) = tmp_Bf; 
clear tmp_Bf Bf_old ind_outlier; 

[u, s, v] = svdsecon(bsxfun(@minus, Bf, mean(Bf, 2)), nb); 
b = u(ind_patch, :)*s; 
f = v'; 
b0 = Ymean(ind_patch) - A_patch(ind_patch,:)*Cmean-b*mean(f, 2);


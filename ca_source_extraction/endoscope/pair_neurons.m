function results = pair_neurons(A1, C1, A2, C2)
%% pair neurons bwtween two demixing results 

%% normalize temporal traces 
C1 = bsxfun(@times, C1, 1./sqrt(sum(C1.^2, 2))); 
C2 = bsxfun(@times, C2, 1./sqrt(sum(C2.^2, 2))); 
A1 = bsxfun(@times, A1, 1./sqrt(sum(A1.^2, 1))); 
A2(A2<0) = 0; 
A2 = bsxfun(@times, A2, 1./sqrt(sum(A2.^2, 1))); 
K = size(A2, 2); 
K1 = size(A1,2); 

results = struct(); 
%% temporal correlation 
C_sim = C2 * C1'; 
[~, results.ind_temporal] = max(C_sim, [], 1); 

%% spatial correaltion 
IND = (A1>1e-5); % spatial mask 
A2norm = sqrt((A2.^2)' * IND) ; 
A2norm(A2norm<1e-5) = inf; 
A_sim = (A2' * A1)./A2norm; 
[~, results.ind_spatial] = max(A_sim, [], 1); 

%% integrated corrlation 
All_sim = A_sim.*C_sim; 
ind1 = bsxfun(@eq, All_sim, max(All_sim, [], 1)); 
ind2 = bsxfun(@eq, All_sim, max(All_sim, [], 2)); 
[val_max, ind_max] = max(ind1.*ind2, [], 1);  % the best match  

ind = sub2ind([K, K1], ind_max, (1:K1)); 
max_spatial = nan(1, K1); 
max_temporal = nan(1, K1); 
max_all = nan(1, K1); 

max_spatial(val_max>0) = A_sim(ind(val_max>0)); 
max_temporal(val_max>0) = C_sim(ind(val_max>0)); 
max_all(val_max>0) = All_sim(ind(val_max>0)); 
ind_max(val_max<=0) = nan; 

results.ind_max = ind_max; 
results.max_spatial = max_spatial; 
results.max_temporal = max_temporal; 
results.max_all = max_all; 
%% results 
function C_ = decorrTemporal(obj, wd)
if ~exist('wd', 'var') || isempty(wd)
    wd = 1; 
end
%% find neighbors of each neurons
A_ = obj.A;
S_ = full(obj.S);
neuron_sn = GetSn(obj.C_raw);
pvals = obj.P.kernel_pars;
S_ = bsxfun(@times, S_, 1./neuron_sn);
ctr = obj.estCenter(); 
temp = sqrt(bsxfun(@minus, ctr(:, 1), ctr(:,1)').^2 + bsxfun(@minus, ctr(:,2), ctr(:,2)').^2); 
ind_neigh = (temp<obj.options.gSiz); 

[K, T] = size(S_);
C_ = obj.C;

num_per_row = 100;
for m=1:K
    fprintf('|');
    if mod(m, num_per_row)==0
        fprintf('\n');
    end
end
fprintf('\n');


%% get kernels
model = obj.options.deconv_options.type;

%% remove small spikes
Tk = 500; 
for k=1:K
    s = S_(k, :);
    tmpS = S_(ind_neigh(k, :), :);
    if isempty(tmpS)
        continue;
    else
        ind = (s~=max(tmpS, [], 1));
    end
    
    if wd>1
        ind = (conv(double(ind), ones(1, wd), 'same') > 0); 
    end
    s(ind) = 0;
    
    fprintf('.');
    
    % generate calcium trace
    if strcmpi(model, 'ar1')
        ht = (pvals(k)).^(0:(Tk-1));
    elseif strcmpi(model, 'ar2')
        ht = exp2kernel(ar2exp(pvals(k,:)), Tk);
    end
    c = conv(double(full(s')), ht/ht(1))*neuron_sn(k);
    C_(k, :) = c(1:T);
    if mod(k,num_per_row)==0
        fprintf('\n');
    end
end
fprintf('\n');


obj.C = C_;
end

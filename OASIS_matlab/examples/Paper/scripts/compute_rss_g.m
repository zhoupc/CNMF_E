function rss_vec = compute_rss_g(g, y, Aset, lam)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   Aset: npools*4 matrix, previous active sets
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, curret value of sparsity penalty parameter lambda.

%% outputs
%   rss_vec: vector with the same elements as g. 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016

%% initialization
len_active_set = size(Aset, 1);  %number of active sets
y = reshape(y,[],1);    % fluorescence data
maxl = max(Aset(:, 4));   % maximum ISI
c = zeros(size(y));     % the optimal denoised trace

%% function for computing the optimal RSS given fixed AR coefficient g and the active set
    function rss = rss_g(g)
        yp = y - lam*(1-g);     % include the penalty term
        h = exp(log(g)*(0:maxl)');   % response kernel
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        for ii=1:len_active_set
            li = Aset(ii, 4);
            ti = Aset(ii, 3);
            idx = ti:(ti+li-1);
            tmp_v = yp(idx)' * h(1:li) / hh(li);
            c(idx) = tmp_v*h(1:li);
        end
        res = y-c;
        rss = res'*res;     % residual sum of squares
    end

rss_vec = zeros(size(g));
for m=1:length(g)
    rss_vec(m) = rss_g(g(m));
end
end

























function vmax = max_ht(pars1, pars2)
%% Get the maximum response value following one calcium transient
%% inputs:
% pars1: if nargin<2, pars1 is the AR coefficients; else pars1 is the
% decay time constant
% pars2: the rising time constant
%% outputs:
%   ht:     nmax*1 kernel, convolution kernel

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
if nargin<2
    if length(pars1) ==1
        % AR1 model
        vmax = 1;
        return;
    else
        % AR2 model
        taus = ar2exp(pars1);
    end
else
    taus = [pars1, pars2];
end

t = (1:ceil(taus(1)*2))';
d = exp(-1./taus(1));
r = exp(-1./taus(2));
ht = (exp(log(d)*t) - exp(log(r)*t) ) / (d - r);
vmax = max(ht);
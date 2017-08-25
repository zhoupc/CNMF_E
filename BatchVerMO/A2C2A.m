function ACS_temp=A2C2A(Ysignal, A, options)
% ACS=A2C2A(ACS, File, A, options)
% General Description:
%      This function is designed to extract C in background-subtracted signal
%      from A(not Amask) in Ysignal.
%      The method relies on CNMF-E's conference script.
%      For each successful ai extraction, peel A*C off Ysignal.
%      Iterate for all the neurons indicated by A.
% Input: 
% Ysignal(background-subtracted), with dimenstion of [d1*d2,T];
% A is also used for center calculation(for subtraction previous signal).
% options: struct data of paramters/options
    %    d1:     number of rows
    %    d2:     number of columns
    %    gSiz:   maximum size of a neuron
    %    nb:     left-over from the previous version, number of knots in
    %            modeling background,=1.
    %    min_corr: minimum threshold of correlation for segementing neurons
    %    deconv_options
    %    deconv_flag    
% Output:
% ACS_temp is a structure with three fields: A,C,and C's
%       standard deviation(STD). This structure is used for merging in later steps in cnmf_e-BatchVer. 

% Modified from "greedyROI_endoscope"
% Main dependences are "extract_a", "extract_c" in CNMF-E BatchVer, and "com" for neuron center by Pengcheng Zhou.

% Shijie Gu, techel@live.cn

%% parameters and preparations

d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSiz = options.gSiz;    % average size of neurons

Acenter = round(com(A, d1, d2)); %nr x 2 matrix, with the center of mass coordinates

Ysignal(isnan(Ysignal)) = 0;     % remove nan values
Ysignal = double(Ysignal);
T = size(Ysignal, 2);

deconv_options_0= options.deconv_options;
deconv_flag = options.deconv_flag;

K = size(A,2);
Cin = zeros(K, T);      % temporal components
Cin_raw = zeros(K, T);  % temporal components
STD = zeros(1,K);       % standard deviation

%% start initialization
for k = 1:K
    r=Acenter(k,1);
    c=Acenter(k,2);
    
    % Use center point data to roughly check whether this is a good starting point
%    ind_p=sub2ind([d1 d2], r, c);
%     y0 = Ysignal(ind_p, :);    
%     y0_std = std(diff(y0));
%     if max(diff(y0))< y0_std % signal is weak
%         Ain(:,k)=A(:,k);       % If it is a "poor" ci, no problem, low STD will let not it contribute much to finalA.
%         Cin(k,:)=y0;           % Since poor ci will get poor A that has bad shapes. Use normal A to fill in the place.
%         STD(k)=std(y0);
%         continue;
%     end
    
    % extract ci    
    [ci_raw,ind_success_ci] = extract_c(Ysignal,[],A(:,k));
    Cin_raw(k,:)=ci_raw;
    
    if ~ind_success_ci          % If it is a "poor" ci, no need to subtract it before next neuron.
        Cin(k,:)=ci_raw;
        STD(k)=std(ci_raw);
        continue;
    end
    
    if and(ind_success_ci,deconv_flag)
        try
            % deconv the temporal trace
            [ci, ~, ~] = deconvolveCa(ci_raw, deconv_options_0);
            % save this initialization
            Cin(k, :) = ci;            
            STD(k)=std(ci);
            ci=ci';
        catch
            Cin(k, :)=ci_raw;
            STD(k)=std(ci_raw);
            ci=ci_raw;
        end      
    end
    
    % select its neighbours to subtract some data, the box size is
    %[2*gSiz+1, 2*gSiz+1]
    % if the c is poor
    ci_std = std(diff(ci));
    if max(diff(ci))< ci_std % signal is weak
        continue
    end
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    HY_box = Ysignal(ind_nhood, :);      
    Ysignal(ind_nhood, :) = HY_box - A(ind_nhood,k)*ci_raw;  % update data
end
ACS_temp=struct('Cin',[],'Cin_raw',[],'STD',[]);
ACS_temp.Cin = Cin;
ACS_temp.Cin_raw=Cin_raw;
ACS_temp.STD = STD;

end
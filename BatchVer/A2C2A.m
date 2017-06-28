function ACS=A2C2A(ACS, File, A, options)
% ACS=A2C2A(ACS, File, A, options)
% General Description:
%      This function is designed to go in de-noised signal
%      from A's Amask to C in File and its paired A.
%      For each neuron, from Amask, get C from mean(see help extract_c), then regress for A.
%      For successful ai extraction, peel A*C off Ysignal.
%      Iterate for all the neurons indicated by Amask.
% Input: 
% ACS is a structure for each file in the folder with fields: A,C,whose
%       standard deviation,STD for all the A's concatnated from each file in
%       the sampling stage. This structure is used for merging. Input here
%       would be an empty but pre-structured structure.
% File is a structure that contains at least one field: Ysignal(denoised-Bgsubtracted)
% A is used for Amask and center calculation(for rough checking and trimming A).
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
% Filling ACS.
% Modified from "greedyROI_endoscope"
% Main dependences are "extract_a", "extract_c", and "com" for neuron center by Pengcheng Zhou.

% Shijie Gu, techel@live.cn

%% parameters and preparations

d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSiz = options.gSiz;    % average size of neurons

Acenter = round(com(A, d1, d2)); %nr x 2 matrix, with the center of mass coordinates
Amask=A>0;

Ysignal=File.Ysignal;
Ysignal(isnan(Ysignal)) = 0;    % remove nan values
Ysignal = double(Ysignal);
T = size(Ysignal, 2);

deconv_options_0= options.deconv_options;
deconv_flag = options.deconv_flag;

K = size(Amask,2);
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
STD = zeros(1,K);       % standard deviation

%% start initialization
for k = 1:K;
    r=Acenter(k,1);
    c=Acenter(k,2);
    
    % Use center point data to roughly check whether this is a good starting point
    ind_p=sub2ind([d1 d2], r, c);
    y0 = Ysignal(ind_p, :);    
    y0_std = std(diff(y0));
    if max(diff(y0))< 3*y0_std % signal is weak
        Ain(:,k)=A(:,k);       % If it is a "poor" ci, no problem, low STD will let not it contribute much to finalA.
        Cin(k,:)=y0;           % Since poor ci will get poor A that has bad shapes. Use normal A to fill in the place.
        STD(k)=std(y0);
        continue;
    end
    
    % extract ci    
    [ci_raw,ind_success_ci] = extract_c(Ysignal,Amask(:,k),[]);    
    
    if ~ind_success_ci   % If it is a "poor" ci, no problem, low STD will let not it contribute much to finalA.
        Ain(:,k)=A(:,k); % Since poor ci will get poor A that has bad shapes. Use normal A to fill in the place.
        Cin(k,:)=ci_raw;
        STD(k)=std(ci_raw);
        continue;
    end
    
    if and(ind_success_ci,deconv_flag)
        try
            % deconv the temporal trace
            [ci, ~, ~] = deconvolveCa(ci_raw, deconv_options_0, 'sn', 1);  % sn is 1 if Ysignal has no noise.
            % save this initialization
            Cin(k, :) = ci;            
            STD(k)=std(ci);
            ci=ci';
        catch
            Cin(k, :)=ci_raw;
            ci=ci_raw';
            STD(k)=std(ci_raw);                        
        end      
    end
    % extract ai        
    % select its neighbours for estimation of ai, the box size is
    %[2*gSiz+1, 2*gSiz+1]
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);               % size of the box
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    HY_box = Ysignal(ind_nhood, :);      
    Amask_box=Amask(ind_nhood,k);
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % index of the center

    sz = [nr, nc];
    [ai, ind_success_ai] = extract_a(ci, [], HY_box, Amask_box, ind_ctr, sz, 1, []); 
    if ind_success_ai==false
        Ain(:,k)=A(:,k);
    else
        Ain(ind_nhood,k)=ai;
        Ysignal(ind_nhood, :) = HY_box - ai*ci;  % update data
    end  
end
    ACS.Ain = [ACS.Ain Ain];
    ACS.Cin = [ACS.Cin; Cin];
    ACS.STD=[ACS.STD STD];
end
function [gSig, gSiz, ring_radius, min_corr, min_pnr] = set_parameters(obj)
%% choose parameters used in CNMF-E

%% Author: Pengcheng Zhou, Columbia University, 2017
%zhoupc1988@gmail.com

gSiz = obj.options.gSiz;
gSig = obj.options.gSig;
ring_radius = obj.options.ring_radius;
min_corr = obj.options.min_corr;
min_pnr = obj.options.min_pnr;

%% choose gSiz, which is the neuron diameter
fprintf('\n*** Determine the neuron diameter (gSiz) and the gaussian width (gSig) of the filtering kernel.\n');
fprintf('\t-------------------------- GUIDE --------------------------\n');
fprintf('\tgSiz is usually slightly larger than a neuron;\n');
fprintf('\tgSig is usually selected as 1/4 of gSiz when the data is 1p\n');
fprintf('\tgSig is usually selected as 1 when the data is 2p.\n');
fprintf('\tIntegers are preferred for both values, but not necessary.\n');
fprintf('\tOdd number is preferred for gSiz. \n');
fprintf('\tWhen you want to make a change, an image will be shown.\n');
fprintf('\tYou can zoon in to see the size of a typical neuron.\n');
fprintf('\t--------------------------  END  --------------------------\n\n');

fprintf('You current selection of parameters (gSiz, gSig) are (%.1f, %.1f).\n', gSiz, gSig);
temp = input('Do you want to make a change? (y/n)    ', 's');
if strcmpi(temp, 'n')
    fprintf('Your values for (gSiz, gSig) will stay the same\n');
else
    Y = get_patch_data(obj.P.mat_data, [], [1,100]);
    [d1,d2, T] = size(Y);
    Y = reshape(Y, d1*d2, T);
    fprintf('Wait for a few second ....\n');
    if obj.options.center_psf
    [u, v] = nnmf(double(Y),1);
        img = reshape(max(double(Y)-u*v, [], 2), d1,d2);
    else
            img = reshape(max(Y,[],2), d1,d2);
    end
    fig_1 = figure;
    imagesc(img); axis equal;
    while true
        temp = input('type new values for (gSiz, gSig):  ', 's');
        try
            new_values = str2num(temp); %#ok<*ST2NM>
            gSiz = new_values(1);
            gSig = new_values(2);
            obj.options.gSig = gSig;
            obj.options.gSiz = gSiz;
            fprintf('Good! You current selection of parameters (gSiz, gSig) are (%.1f, %.1f).\n', gSiz, gSig);
            break;
        catch
            warning('values are bad. try again');
            continue;
        end
    end
    try
        close(fig_1);
    end
end

%%
%% choose the radius of the ring
if strcmpi('obj.options.background_model', 'ring')
    fprintf('\n*** Determine the ring radius for estimating background fluctuations.\n');
    fprintf('\t-------------------------- GUIDE --------------------------\n');
    fprintf('\tThe ring radius is selected to be larger than neuron diameters. \n');
    fprintf('\tThus pixels in the center only shares the common background fluctuations with pixels on the ring.\n');
    fprintf('\tA usual ratio between ring radius and neuron diameter is between 1 to 3. 1.5 is the default value.\n');
    fprintf('\tSmall values have a change of being smaller than neuron size.\n');
    fprintf('\tLarge values increase the computation cost.\n');
    fprintf('\t--------------------------  END  --------------------------\n\n');
    
    fprintf('You current selection of ring_radius is %d pixels.\n', ring_radius);
    temp = input('Do you want to make a change? (y/n)    ', 's');
    if strcmpi(temp, 'n')
        fprintf('Your values for ring_radius will stay the same\n');
    else
        while true
            temp = input('type new values for ring_radius:  ', 's');
            try
                ring_radius = str2double(temp);
                obj.options.ring_radius = ring_radius;
                fprintf('Good! You current selection of parameters ring_radius is %2d.\n', ring_radius);
                break;
            catch
                warning('values are bad. try again');
                continue;
            end
        end
    end
end
%% choose parameters for min_corr and min_pnr
fprintf('\n*** Determine min_corr and min_pnr for screening seed pixels.\n');
fprintf('\t-------------------------- GUIDE --------------------------\n');
fprintf('\tA seed pixel is used for initializing one neuron.\n');
fprintf('\tSeed pixel usually locates at the center of each neuron.\n');
fprintf('\tSeed pixel is determined using two values: min_corr and min_pnr.\n');
fprintf('\tmin_pnr stands for minimum local correlation. \n');
fprintf('\tPNR stands for minimum peak-to-noise ratio.\n');
fprintf('\tWe want to pick thresholds that all neuron centers are above the thresholds,\n');
fprintf('\tand include the least number of non-ROI pixels. \n');
fprintf('\tThere is a trade-of between missing neurons and acquiring more false positives.\n');
fprintf('\tWhen it is hard to determine, we tend to choose missing neurons, i.e., larger threshlds.\n');
fprintf('\tThose missed neurons can be picked up from the residual video later.\n');
fprintf('\tTip: use Data Cursor to help you pick thresholds.\n');
fprintf('\t--------------------------  END  --------------------------\n\n');

fprintf('You current selection of parameters (min_corr, min_pnr) are (%.3f, %.3f).\n', min_corr, min_pnr);
temp = input('Do you want to make a change? (y/n)    ', 's');
if strcmpi(temp, 'n')
    fprintf('Your values for (min_corr, min_pnr) will stay the same\n');
else
    
    if isempty(obj.Cn) || isempty(obj.PNR)
        fprintf('computing the correlation image and the peak-to-noise ratio image....\n');
        [cn, pnr] = obj.correlation_pnr_parallel([1, 5000]);
        obj.Cn = cn;
        obj.PNR = pnr;
        fprintf('Done\n');
    else
        cn = obj.Cn;
        pnr = obj.PNR;
    end
    
    %find all local maximum as initialization point
    tmp_d = max(1,round(gSiz/4));
    v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
    ind = (v_max==cn.*pnr);
    
    figure('papersize', [12, 3]);
    init_fig;
    subplot(131);
    imagesc(cn, [0, 1]);
    title('local corr. image');
    axis equal off tight;
    subplot(132);
    pnr_vmax = max(pnr(:))*0.8;
    imagesc(pnr, [3, pnr_vmax]);
    axis equal off tight;
    title('PNR image');
    subplot(133);
    imagesc(cn.*pnr, [3, pnr_vmax]);
    hold on;
    
    tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
    [r, c] = find(tmp_ind);
    ax_seeds = plot(c, r, '.m', 'markersize', 10);
    axis equal off tight;
    title('candidate seed pixels');
    ylabel('PNR');
    xlabel('Cn');
    
    while true
        temp = input('type new values for (min_corr, min_pnr):  ', 's');
        try
            new_values = str2num(temp); %#ok<*ST2NM>
            min_corr = new_values(1);
            min_pnr = new_values(2);
            obj.options.min_corr = min_corr;
            obj.options.min_pnr = min_pnr;
            tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
            [r, c] = find(tmp_ind);
            delete(ax_seeds);
            subplot(133); 
            ax_seeds = plot(c, r, '.m', 'markersize', 10);
            
            temp = input('Is this good? (y/n)   ', 's');
            if strcmpi(temp, 'y')
                fprintf('Good! You current selection of parameters (min_corr, min_pnr) are (%2d, %2d).\n', min_corr, min_pnr);
                break;
                
            end
        catch
            warning('values are bad. try again');
            continue;
        end
    end
end

fprintf('\nThe parameters used in CNMF-E are \ngSig:\t\t%d\ngSiz:\t\t%d\nring_radius\t%d\nmin_corr:\t%.3f\nmin_pnr:\t%.3f\n', ...
    gSig, gSiz, ring_radius, min_corr, min_pnr);
end
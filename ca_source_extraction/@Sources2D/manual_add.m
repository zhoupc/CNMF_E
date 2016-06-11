function manual_add(obj, Yres)
%% manually add missing neurons

if ~ismatrix(Yres); Yres = obj.reshape(Yres, 1); end;
Ymax = max(Yres, [], 2);
d1 = obj.options.d1;
d2 = obj.options.d2;
gSiz = obj.options.gSiz;
gSig = obj.options.gSig;

figure;
subplot(221);
obj.image(Ymax);
axis equal off tight;
POS = get(gca, 'position'); 
POS(1) = POS(1)+POS(3)+0.01; 
POS(3) = 0.02; 
S='temp=subplot(221); imagesc(temp.Children.CData, [0, get(gco, ''value'')]); axis equal off tight;';
uic=uicontrol('style','slider','units','normalized','position', ...
    POS,'min',0,'max',max(Ymax(:)),'value',max(Ymax),'callback',S);

A = zeros(size(obj.A));
C = zeros(size(obj.C));
k = 0;  % number of neurons added
psf = fspecial('gaussian', round(gSiz), gSig);
psf = psf-mean(psf(:));
disp('adjust the colorbar and type any key to continue\n'); 
pause; 

while true
    % pick neuron locations
    subplot(221);
    obj.image(Ymax, [0, get(uic, 'value')]);
    axis equal off tight; 
    [c,r] = ginput(1); c=round(c); r=round(r);
    ind_ctr = sub2ind([d1,d2], r, c);
    y0 = squeeze(Yres(ind_ctr, :));
    
    % select neighbours
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    Y_box = Yres(ind_nhood, :);
    HY_box = imfilter(reshape(Y_box, nr, nc,[]), psf, 'replicate');
    HY_box = reshape(HY_box, nr*nc, []);
    
    % run rank-1 NMF
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
    sz = [nr, nc];
    [ai, ci, ind_success] =  extract_ac(HY_box, Y_box, ind_ctr, sz);
    if ~ind_success
        continue;
    end
    subplot(222);
    imagesc(reshape(ai, nr, nc)); axis equal off tight;
    subplot(2,2,3:4);
    plot(ci);
    
    temp = questdlg('save this neuron? ', '', 'No');
    if strcmpi(temp, 'Yes')
        k = k+1;
        A(ind_nhood, k) = ai;
        C(k, :) = ci;
    elseif strcmpi(temp, 'Cancel')
        break;
    end
    Yres(ind_nhood, :) = Y_box - ai*ci;
    Ymax(ind_nhood) = max(Yres(ind_nhood, :), [], 2);
end

if k~=0
    obj.A = [obj.A, A(:, 1:k)];
    obj.C = [obj.C; C(1:k, :)];
end

























function  Brainbow=ColorAllNeurons(A,d1,d2,Picname,outputdir)
    if nargin <5
        outputdir = {};
    end
    if nargin<4
        Picname = {};
    end

%     figure; 
    hold all
    color_palet = 1-[[1 0 0]; [1 .6 0]; [.7 .6 .4]; [.6 .8 .3]; [0 .6 .3]; [0 0 1]; [0 .6 1]; [0 .7 .7]; [.7 0 .7];  [.7 .4 1]]; 
    color_palet = color_palet([1:2:end 2:2:end],:); % scramble slightly
    nColors = color_palet(mod(1:(size(A,2)),size(color_palet,1))+1,:); 
    Ir = diag(nColors(:,1)); 
    Ig = diag(nColors(:,2)); 
    Ib = diag(nColors(:,3));
    
    
    k = size(A,2); 
    sA = sum(A,1);
    A = bsxfun(@rdivide, A, sA)*prctile(sA, 5);
    
    C=ones(k,1);
    
    Brainbow = 1-cat(3, A*Ir*C, A*Ig*C, A*Ib*C); 
    Brainbow = reshape(Brainbow,d1,d2,3);
    image(Brainbow)
    set(gca, 'ydir', 'reverse')
    axis image; axis off; shg
    
    Atemp=reshape(A,d1,d2,k);
    Position=zeros(2,k);
    for i=1:k
        Atemp_=Atemp(:,:,i);
        [row_ind,col_ind] = find(Atemp_>0);
        Position(2,i)=mean(row_ind);
        Position(1,i)=mean(col_ind);
    end    
    text(Position(1,:),Position(2,:),cellstr(num2str((1:k)'))','Color','black')
    
    fig=gcf;
    if numel(Picname)>0
        title(Picname,'interpreter','none')

        fignam=[outputdir Picname,'.png'];

        saveas(gcf,fignam);
    end
    
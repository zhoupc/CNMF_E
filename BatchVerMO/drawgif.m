function drawgif(M3,d1,d2,filename,names)
if nargin<5
    names=1:numel(M3);
end
figure
    for i = 1:numel(M3)
        A=M3{i};
        ColorAllNeurons(A,d1,d2,['Day ' num2str(names(i))],[]);
        h=getframe(1);
        [A,map] = rgb2ind(h.cdata,256);
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end
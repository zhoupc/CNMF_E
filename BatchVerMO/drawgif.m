function drawgif(M3,d1,d2,filename)
    for i = 1:numel(M3)
        A=M3{i};
        ColorAllNeurons(A,d1,d2,num2str(i),[]);
        h=getframe(1);
        [A,map] = rgb2ind(h.cdata,256);
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end
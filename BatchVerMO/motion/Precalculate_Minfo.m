function handles=Precalculate_Minfo(M,handles)
diagr=tril(ones(numel(M),numel(M)),-1);
[row,col]=find(diagr);
handles.allpairs=[col,row];

% Opening look
% axes-creat the overlay picture
handles=Update_Plot(1,handles);
linkaxes([handles.axes1 handles.axes2])

% Opening look

% % slider1
set(handles.slider1,'Enable','off')
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',size(handles.allpairs,1));
set(handles.slider1,'Value',1);
set(handles.slider1,'SliderStep', [1/(size(handles.allpairs,1)-1), 1/(size(handles.allpairs,1)-1)]);
function varargout = ManualShift(varargin)
% MANUALSHIFT MATLAB code for ManualShift.fig
%      MANUALSHIFT, by itself, creates a new MANUALSHIFT or raises the existing
%      singleton*.
%
%      H = MANUALSHIFT returns the handle to a new MANUALSHIFT or the handle to
%      the existing singleton*.
%
%      MANUALSHIFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALSHIFT.M with the given input arguments.
%
%      MANUALSHIFT('Property','Value',...) creates a new MANUALSHIFT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ManualShift_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ManualShift_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ManualShift

% Last Modified by GUIDE v2.5 07-Aug-2017 14:11:47


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ManualShift_OpeningFcn, ...
                   'gui_OutputFcn',  @ManualShift_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ManualShift is made visible.
function ManualShift_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ManualShift (see VARARGIN)

% Parsing inputs
guidata(hObject, handles);

p = inputParser;
addParameter(p,'M',[]);
addParameter(p,'d1',300);
addParameter(p,'d2',400);

parse(p, varargin{:});
handles.M=p.Results.M;
handles.M_intact=p.Results.M; % for some cancel buttons.
handles.d1=p.Results.d1;
handles.d2=p.Results.d2;
handles.pickOrnot=false;

% All pairs of alignment
M=p.Results.M;
diagr=tril(ones(numel(M),numel(M)),-1);
[row,col]=find(diagr);
handles.allpairs=[col,row];


% Opening look
% axes-creat the overlay picture
handles=Update_Plot(1,handles);

linkaxes([handles.axes1 handles.axes2])


% % slider1
set(handles.slider1,'Enable','off')
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',size(handles.allpairs,1));
set(handles.slider1,'Value',1);
set(handles.slider1,'SliderStep', [1/(size(handles.allpairs,1)-1), 1/(size(handles.allpairs,1)-1)]);

% Choose default command line output for ManualShift
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ManualShift wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualShift_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
display('Saving results...')
varargout{1} = handles.output;
varargout{2} = handles.M;
varargout{3} = handles.M_intact;
delete(handles.figure1);

function neuron_list_tmpl_Callback(hObject, eventdata, handles)
% hObject    handle to neuron_list_tmpl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron_list_tmpl as text
%        str2double(get(hObject,'String')) returns contents of neuron_list_tmpl as a double
neuron_list_str=get(hObject,'String');
display(neuron_list_str)
if ~strcmp(neuron_list_str(end),';')
    neuron_list_str(end+1)=';';
end
neuron_list=sscanf(neuron_list_str,'%d-%d;');
handles.neuron_list=reshape(neuron_list,2,[]);  %two rows, first row-template, second row-the toAlign.
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function neuron_list_tmpl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron_list_tmpl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

currentpair_ind=get(hObject,'Value');
currentpair_ind = round(currentpair_ind); %round off this value
set(hObject, 'Value', currentpair_ind);

handles=Update_Plot(currentpair_ind,handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',1);


function ToAlign_Callback(hObject, eventdata, handles)
% hObject    handle to ToAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToAlign as text
%        str2double(get(hObject,'String')) returns contents of ToAlign as a double
ToAlign_text=get(hObject,'String');
ToAlign_num=sscanf(ToAlign_text,'%d-%d');
handles.currentpair(1)=ToAlign_num(1);
handles.currentpair(2)=ToAlign_num(2);
handles=Update_Plot([],handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ToAlign_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ToAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','1-2')


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neuron_list_tmpl_Callback(handles.neuron_list_tmpl,eventdata, handles)
handles=guidata(hObject);
d1=handles.d1;                             d2=handles.d2;
template_num=handles.currentpair(1);       ToAlign_num=handles.currentpair(2);
neuron_templ_ind=handles.neuron_list(1,:); neuron_ToAlign_ind=handles.neuron_list(2,:);

for i = 1:length(neuron_templ_ind)
    template_neuron_A=handles.M{template_num}{template_num}(:,neuron_templ_ind(i));
    ToAlign_neuron_A=handles.M{template_num}{ToAlign_num}(:,neuron_ToAlign_ind(i));
            % template_neuron_center = round(com(template_neuron_A, d1, d2));
            % ToAlign_neuron_center = round(com(ToAlign_neuron_A, d1, d2));
            % displacement_vector=[ToAlign_neuron_center template_neuron_center];
    [D,~] = imregdemons(A2image(ToAlign_neuron_A,d1,d2),A2image(template_neuron_A,d1,d2));
    handles.M{template_num}{ToAlign_num}(:,neuron_ToAlign_ind(i))=...
        reshape(imwarp(reshape(ToAlign_neuron_A,d1,d2),D,'cubic'),[],1);

    template_neuron_A_reverse=handles.M{ToAlign_num}{template_num}(:,neuron_templ_ind(i));
    handles.M{ToAlign_num}{template_num}(:,neuron_templ_ind(i))=...
        reshape(imwarp(reshape(template_neuron_A_reverse,d1,d2),D.*(-1),'cubic'),[],1);
end
display('A Updated')
guidata(hObject,handles);
ToAlign_Callback(handles.ToAlign, eventdata, handles)
% Update handles structure
%guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
template_num=handles.currentpair(1);       ToAlign_num=handles.currentpair(2);
neuron_templ_ind=handles.neuron_list(1,:); neuron_ToAlign_ind=handles.neuron_list(2,:);

for i = 1:length(neuron_templ_ind)

    handles.M{template_num}{ToAlign_num}(:,neuron_ToAlign_ind(i))=...
            handles.M_intact{template_num}{ToAlign_num}(:,neuron_ToAlign_ind(i));
    handles.M{ToAlign_num}{template_num}(:,neuron_templ_ind(i))=...
            handles.M_intact{ToAlign_num}{template_num}(:,neuron_templ_ind(i));
end

% Update handles structure
guidata(hObject, handles);
ToAlign_Callback(handles.ToAlign, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
set(gca,'tag','axes1')

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
set(gca,'tag','axes2')

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% KeyPressed = eventdata.Key;
% switch KeyPressed
%     case 'rightarrow'
%         axes(handles.axes1)
%         imshow(handles.B)
%     case 'leftarrow'
%         axes(handles.axes1)
%         imshow(handles.A)
%     case 'downarrow'
%         axes(handles.axes1)
%         imshowpair(handles.A,handles.B,'falsecolor','Scaling','independent');
%     case 'uparrow'
%         set(hObject,'CurrentObject',handles.slider1)
% end


% --- Executes on key press with focus on slider1 and none of its controls.
function slider1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
uiresume(hObject);
d1=handles.d1;
d2=handles.d2;
KeyPressed = eventdata.Key;
switch KeyPressed
    case 'space'
        slider_enable=get(handles.slider1,'Enable');
        if strcmp(slider_enable,'on')
            set(handles.slider1, 'Enable', 'off');
        elseif strcmp(slider_enable,'off')
            set(handles.slider1, 'Enable', 'on');
            uicontrol(handles.slider1)
        end        
    case 'rightarrow'
        handles.mode='toAlign';
        axes(gca)
        current_axes=gca;
        if strcmp(current_axes.Tag,'axes1')
            imshow(handles.B)
            CalText_AddText(handles.B_text);
            handles.mode='toAlign';
            set(gca,'tag','axes1')
            set(handles.text6,'String','Registered New A To Template');
        elseif strcmp(current_axes.Tag,'axes2')
            imshow(handles.D)
            CalText_AddText(handles.D_text);
            handles.mode='toAlign';
            set(gca,'tag','axes2');
            set(handles.text7,'String','Original A to Align');
        end
    case 'leftarrow'
        handles.mode='template';
        axes(gca)
        current_axes=gca;
        if strcmp(current_axes.Tag,'axes1')
            imshow(handles.A); CalText_AddText(handles.A_text); 
            set(gca,'tag','axes1'); set(handles.text6,'String','Template');
        elseif strcmp(current_axes.Tag,'axes2')
            imshow(handles.A); CalText_AddText(handles.A_text);
            set(gca,'tag','axes2'); set(handles.text7,'String','Template'); 
            handles.figure1.CurrentAxes=handles.axes2;
        end

    case 'downarrow'
        handles.mode='overlapping';
        axes(gca)
        current_axes=gca;
        if strcmp(current_axes.Tag,'axes1')
            imshow(handles.C);
            CalText_AddText(handles.A_text); CalText_AddText(handles.B_text);
            %imshowpair(handles.A_b,handles.B_b,'falsecolor','Scaling','independent');
            set(gca,'tag','axes1')
            set(handles.text6,'String','Registered(Green) and Template(Red)'); 
        elseif strcmp(current_axes.Tag,'axes2')
            imshow(handles.E);
            CalText_AddText(handles.A_text); CalText_AddText(handles.D_text);
            %imshowpair(handles.A_b,handles.D_b,'falsecolor','Scaling','independent');
            set(gca,'tag','axes2');
            set(handles.text7,'String','Un-Registered(Green) and Template(Red)');
            handles.figure1.CurrentAxes=handles.axes2;
        end        
    case 'uparrow'
        set(hObject,'CurrentObject',handles.slider1)
end
guidata(hObject, handles);
uiwait(hObject);

function handles=Update_Plot(currentpair_ind,handles)
handles.mode='overlapping';
d1=handles.d1;  d2=handles.d2;
if ~isempty(currentpair_ind)
    handles.currentpair=handles.allpairs(currentpair_ind,:);
end
template_num=handles.currentpair(1);
ToAlign_num=handles.currentpair(2);
set(handles.ToAlign,'String',[num2str(template_num) '-' num2str(ToAlign_num)]);
set(handles.text10,'String',[]);
[~,currentpair_ind,~] = intersect(handles.allpairs,handles.currentpair,'rows');
if ~isempty(currentpair_ind)
    set(handles.slider1,'Value',currentpair_ind)
end
set(handles.neuron_list_tmpl,'String',[]);

M=handles.M;
A=A2image(M{template_num}{template_num},d1,d2,false,'magenta');  handles.A=A;
A_b=A2image(M{template_num}{template_num},d1,d2,false); % black and white for imshowpair  
handles.A_b=A_b;
A_text=CalText_AddText(false,M{template_num}{template_num},d1,d2);                           handles.A_text=A_text;

B=A2image(M{template_num}{ToAlign_num},d1,d2,false,'green');  handles.B=B;
B_b=A2image(M{template_num}{ToAlign_num},d1,d2,false);  B_b = imhistmatch(B_b,A_b); handles.B_b=B_b;
B_text=CalText_AddText(false,M{template_num}{ToAlign_num},d1,d2);                        handles.B_text=B_text;

D=A2image(M{ToAlign_num}{ToAlign_num},d1,d2,false,'green');   handles.D=D;
D_b=A2image(M{ToAlign_num}{ToAlign_num},d1,d2,false);   D_b = imhistmatch(D_b,A_b); handles.D_b=D_b;
D_text=CalText_AddText(false,M{ToAlign_num}{ToAlign_num},d1,d2);                        handles.D_text=D_text;

axes(handles.axes1)
% imshowpair(A_b,B_b,'falsecolor','Scaling','independent')
C=imfuse(A_b,B_b,'falsecolor','Scaling','independent');       handles.C=C;
imshow(C);
CalText_AddText(A_text);
CalText_AddText(B_text);
set(gca,'tag','axes1')
axis tight

axes(handles.axes2)
%imshowpair(A_b,D_b,'falsecolor','Scaling','independent')
E=imfuse(A_b,D_b,'falsecolor','Scaling','independent');       handles.E=E;
imshow(E);
CalText_AddText(A_text);
CalText_AddText(D_text);
set(gca,'tag','axes2')
axis tight


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiresume(hObject);
% uicontrol(hObject)
% uiwait(hObject);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(hObject);
%uicontrol(handles.figure1)

    cursorPoint = get(gca, 'CurrentPoint');
    current_axes=gca;
    x=round(cursorPoint(1,1,1)); y=round(cursorPoint(1,2,1));
    display(x)
    display(y)
    d1=handles.d1;                      d2=handles.d2;              M=handles.M;
    if and(x<=d2,y<=d1)
        mouse_location = sub2ind([d1 d2], y, x);
        template_num=handles.currentpair(1);ToAlign_num=handles.currentpair(2);
        if strcmp(handles.mode,'template')
            A_=M{template_num}{template_num};
        elseif strcmp(handles.mode,'toAlign')
            if strcmp(current_axes.Tag,'axes1')
                A_=M{template_num}{ToAlign_num};
            elseif strcmp(current_axes.Tag,'axes2')
                A_=M{ToAlign_num}{ToAlign_num};
            end
        else
            disp('You are in "overlapping" mode.')
        end
        try
            neuron_ind=find(A_(mouse_location,:)>0);
            set(handles.text10,'String',['neuron #' num2str(neuron_ind)])
        catch
            disp('Neuron picking is not supported in overlapping mode.')
        end
        %guidata(hObject, handles);
    end
%end
uiwait(hObject);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end

function neuron_list_tmpl_KeyPressFcn(hObject, eventdata, handles)
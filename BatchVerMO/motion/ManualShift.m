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

% Last Modified by GUIDE v2.5 03-Aug-2017 18:14:49


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
p = inputParser;
addParameter(p,'M',[]);
parse(p, varargin{:});
handles.M=p.Results.M; 

% All pairs of alignment
M=p.Results.M;
diagr=tril(ones(numel(M),numel(M)),-1);
[row,col]=find(diagr);
handles.allpairs=[col,row];
handles.currentpair=[1,2];

template_num=handles.currentpair(1);
ToAlign_num=handles.currentpair(2);

% Creat the overlay picture
A=A2image(M{template_num}{template_num},300,400);
B=A2image(M{template_num}{ToAlign_num},300,400);
B = imhistmatch(B,A);
imshowpair(A,B,'falsecolor','Scaling','independent');

axes(handles.axes2)
imshow(A)

axes(handles.axes3)
imshow(B)
linkaxes([axes1 axes2 axes3])


% Choose default command line output for ManualShift
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ManualShift wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualShift_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function TemplateChoice_Callback(hObject, eventdata, handles)
% hObject    handle to TemplateChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TemplateChoice as text
%        str2double(get(hObject,'String')) returns contents of TemplateChoice as a double
cuurent_template=get(hObject,'String');


% --- Executes during object creation, after setting all properties.
function TemplateChoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TemplateChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function neuron_list_tmpl_Callback(hObject, eventdata, handles)
% hObject    handle to neuron_list_tmpl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron_list_tmpl as text
%        str2double(get(hObject,'String')) returns contents of neuron_list_tmpl as a double


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



function neuron_list_2align_Callback(hObject, eventdata, handles)
% hObject    handle to neuron_list_2align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of neuron_list_2align as text
%        str2double(get(hObject,'String')) returns contents of neuron_list_2align as a double


% --- Executes during object creation, after setting all properties.
function neuron_list_2align_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron_list_2align (see GCBO)
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
handles.currentpair=handles.allpairs(currentpair_ind,:);
template_num=handles.currentpair(1);
ToAlign_num=handles.currentpair(2);

axes(handles.axes1)
A=A2image(M{template_num}{template_num},300,400);
B=A2image(M{template_num}{ToAlign_num},300,400);
imshowpair(A,B,'falsecolor','Scaling','independent');

axes(handles.axes2)
imshow(A)

axes(handles.axes3)
imshow(B)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',1);
set(hObject,'Max',size(handles.allpairs,1));
handles.currentpair=[1,2];


function ToAlign_Callback(hObject, eventdata, handles)
% hObject    handle to ToAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ToAlign as text
%        str2double(get(hObject,'String')) returns contents of ToAlign as a double


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

%------------suppose there is handles.neuron_templ_ind and
%handles.neuron_ToAlign_ind

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d1=handles.d1;
d2=handles.d2;
template_num=handles.currentpair(1);
ToAlign_num=handles.currentpair(2);
template_neuron_A=handles.M{template_num}{template_num}(:,handles.neuron_templ_ind);
ToAlign_neuron_A=handles.M{template_num}{ToAlign_num}(:,handles.neuron_ToAlign_ind);
% template_neuron_center = round(com(template_neuron_A, d1, d2));
% ToAlign_neuron_center = round(com(ToAlign_neuron_A, d1, d2));
% displacement_vector=[ToAlign_neuron_center template_neuron_center];
[D,~] = imregdemons(A2image(ToAlign_neuron_A),A2image(template_neuron_A));

handles.M{template_num}{ToAlign_num}(:,handles.neuron_ToAlign_ind)=...
    reshape(imwarp(reshape(ToAlign_neuron_A,d1,d2),D,'cubic'),[],1);

template_neuron_A_reverse=handles.M{ToAlign_num}{template_num}(:,handles.neuron_templ_ind);
handles.M{ToAlign_num}{template_num}(:,handles.neuron_templ_ind)=...
    reshape(imwarp(reshape(template_neuron_A_reverse,d1,d2),D.*(-1),'cubic'),[],1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2

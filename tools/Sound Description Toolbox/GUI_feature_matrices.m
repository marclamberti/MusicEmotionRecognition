function varargout = GUI_feature_matrices(varargin)
% GUI_FEATURE_MATRICES M-file for GUI_feature_matrices.fig
%      GUI_FEATURE_MATRICES, by itself, creates a new GUI_FEATURE_MATRICES or raises the existing
%      singleton*.
%
%      H = GUI_FEATURE_MATRICES returns the handle to a new GUI_FEATURE_MATRICES or the handle to
%      the existing singleton*.
%
%      GUI_FEATURE_MATRICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FEATURE_MATRICES.M with the given input arguments.
%
%      GUI_FEATURE_MATRICES('Property','Value',...) creates a new GUI_FEATURE_MATRICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_feature_matrices_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_feature_matrices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help GUI_feature_matrices

% Last Modified by GUIDE v2.5 30-Jun-2005 20:26:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_feature_matrices_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_feature_matrices_OutputFcn, ...
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


% --- Executes just before GUI_feature_matrices is made visible.
function GUI_feature_matrices_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_feature_matrices (see VARARGIN)

% Choose default command line output for GUI_feature_matrices
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_feature_matrices wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_feature_matrices_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
GUI_main;
delete(handles.figure1);
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
format('long');
% Kaleitai to GUI
[m,p] = uigetfile('*.wav','Load one or more wav files (*.wav)','MultiSelect', 'on');
if p == 0  V = []; W = []; H = []; M = []; return; end
a=char(m);
[b1 b2]=size(m);
length(m);
M = [];

%disp('Ypologismos Feature Matrix gia:');
if iscell(m) 
    h = waitbar(0,'Calculating feature matrices, please wait...');
    for i=1:length(m)
        drawnow
        file = [char(p) char(m{i})];
        drawnow
        waitbar((i-1)/length(m));
        drawnow
        save_feature_matrix(file);
        drawnow
        
        cellFile = cellstr(['Completed feature matrix for: ' char(m{i})]);
        M = [M;cellFile];
        waitbar(i/length(m));
        set(handles.listbox1,'String',M);
        drawnow
    end
    close(h);
else
[s1,s2] = size(m);
   for i=1:s1
       a=strcat(m(i,:));
       file = [char(p) a];
       drawnow;
       save_feature_matrix(file);
       
       cellFile = cellstr(['Completed feature matrix for: ' a]);
       M = [M;cellFile];
       set(handles.listbox1,'String',M);
   end
end
msgbox('Completed Creating Feature Matrices!','Done','help');

% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



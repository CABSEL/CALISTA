function varargout = CALISTA_LI_GUI(varargin)
% CALISTA_LI_GUI MATLAB code for CALISTA_LI_GUI.fig
%      CALISTA_LI_GUI, by itself, creates a new CALISTA_LI_GUI or raises the existing
%      singleton*.
%
%      H = CALISTA_LI_GUI returns the handle to a new CALISTA_LI_GUI or the handle to
%      the existing singleton*.
%
%      CALISTA_LI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALISTA_LI_GUI.M with the given input arguments.
%
%      CALISTA_LI_GUI('Property','Value',...) creates a new CALISTA_LI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CALISTA_LI_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CALISTA_LI_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CALISTA_LI_GUI

% Last Modified by GUIDE v2.5 02-Mar-2018 18:03:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CALISTA_LI_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CALISTA_LI_GUI_OutputFcn, ...
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


% --- Executes just before CALISTA_LI_GUI is made visible.
function CALISTA_LI_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CALISTA_LI_GUI (see VARARGIN)

% Choose default command line output for CALISTA_LI_GUI
handles.output = hObject;
% set(handles.axes1, 'units', 'normalized', 'position', [0.05 0.15 0.9 0.8])

% Update handles structure
guidata(hObject, handles);
%  set(handles.figure1,'Units','Pixels','Position',get(0,'ScreenSize'))
INPUTS=evalin('base','INPUTS');

% Upload and pre-process data
[DATA,INPUTS] = import_data(INPUTS);
%% *** 2-SINGLE-CELL CLUSTERING ***
%
% Type 'help CALISTA_clustering_main' for more information.
[Results, DATA, INPUTS]=CALISTA_clustering_main(DATA,INPUTS);
assignin('base','exit_GUI',0);

% Upload DATA and Results to/from Workspace

assignin('base','DATA',DATA);
assignin('base','INPUTS',INPUTS);
assignin('base','Results',Results);
        
%  Cluster removal (if desired)

cluster_cut=input('Press 1 if you want to remove cell cluster(s). Press 0 otherwise: ');

if cluster_cut==1
    
    which_cut=input('Enter cluster index for removal (e.g 1 or [5 3]): ');
   cells_2_cut2=find(ismember(Results.final_groups,which_cut)==1);
   cells_2_cut2=DATA.cut_sort.abs_indices(cells_2_cut2);
   csvwrite('Cells_2_remove',cells_2_cut2)   
   RESULTS.cluster_cut = 1;
   fprintf('Cell''s indices to remove are saved in "Cells_2_remove.csv". Please run MAIN script again \n')
   return
    
else
    
    Proceed=input(' Press 1 if you want to perform additional analysis (e.g. lineage inference, cell ordering) , 0 otherwise: ');
    if ~Proceed
        assignin('base','exit_GUI',1);
        return
    else
        
        % Calculate cluster distances
        % Plot initial clusters and connected edges
        % arrayfun(@cla,findall(0,'type','axes')) %refresh
        set(handles.edgestag,'Visible','on')
        set(handles.checkbox2,'Value',0)
        
        
        % Create a temporarly Results_GUI structure
        Results_GUI=Results;
        DATA=evalin('base','DATA');
        INPUTS=evalin('base','INPUTS');
        handles_graphtextwarning=handles.graphtextwarning;
        set(handles_graphtextwarning,'String','The graph is connected.','ForegroundColor','b')
        handles_pushbutton1=handles.pushbutton1;
        handles_axes1=handles.axes1;
        handles_uipanel1=handles.uipanel1;
        axes(handles_axes1)
        cla(handles_axes1);
        [Results_GUI]=cluster_distance_GUI(DATA,INPUTS,Results_GUI,handles_graphtextwarning,handles_axes1,handles_uipanel1);
        [Results_GUI]=cluster_distance_GUI_update(Results_GUI,handles_axes1);
        % Update the workspace
        assignin ('base','handles_pushbutton1',handles_pushbutton1)
        assignin ('base','handles_graphtextwarning',handles_graphtextwarning)
        assignin ('base','handles_axes1',handles_axes1)
        assignin ('base','handles_uipanel1',handles_uipanel1)
        assignin ('base','Results_GUI',Results_GUI)
    end

end




% UIWAIT makes CALISTA_LI_GUI wait for user response (see UIRESUME)
% uiwait(handles.axes1);


% --- Outputs from this function are returned to the command line.
function varargout = CALISTA_LI_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
exit_GUI=evalin('base','exit_GUI');
if exit_GUI
    delete(hObject)
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% % --- Executes on selection change in listbox1.
% function listbox1_Callback(hObject, eventdata, handles)
% % hObject    handle to listbox1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from listbox1
% 
% 
% % --- Executes during object creation, after setting all properties.
% function listbox1_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to listbox1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% Results=evalin('base','Results');
% set(hObject,'Max',size(Results.list_box_text,1),'Min',2);
% 
Results_GUI=evalin('base','Results_GUI');
DATA=evalin('base','DATA');
INPUTS=evalin('base','INPUTS');
hh=Results_GUI.TRANSITION.final_graph;
h=hh;
nodes_connection3=[hh.Edges.EndNodes(:,1) hh.Edges.EndNodes(:,2)];
Results_GUI.TRANSITION.nodes_connection=nodes_connection3;
Results=Results_GUI;



assignin ('base','Results_GUI',Results_GUI)
Results= after_LI_GUI(DATA,Results);
assignin ('base','Results',Results)
fprintf('\nFinal lineage progression graph saved.\n')


%% *** 4-DETERMINATION OF TRANSITION GENES ***
%
% Type 'help CALISTA_transition_genes_main' for more information.

[Results]=CALISTA_transition_genes_main(DATA,INPUTS,Results);
assignin ('base','Results',Results)
%% *** 5-PSEUDOTEMPORAL ORDERING OF CELLS ***
%
% Type 'help CALISTA_ordering_main' for more information.

[Results]=CALISTA_ordering_main_2(DATA,INPUTS,Results);
assignin ('base','Results',Results)
%% *** 6-LANDSCAPE PLOTTING ***
%
% Type 'help CALISTA_landscape_plotting_main' for more information.

CALISTA_landscape_plotting_main(INPUTS,Results);


%% *** 7-PATH ANALYSIS ***
%
% Type 'help CALISTA_path_main' for more information.
fprintf('Press 1 if you want to select transition paths and perform additional analysis, 0 otherwise: ')
Proceed2=input('');

if ~Proceed2
    assignin('base','exit_GUI',1);
    return
else
    fprintf('Press 1 if you want to manually upload the selected gene list, 0 otherwise: ')
    Proceed3=input('');
    if ~Proceed3
        Results=CALISTA_path_main(DATA,INPUTS,Results);
        assignin ('base','Results',Results)
    else
        [FileName,PathName,FilterIndex] = uigetfile('*.*');
        filename=strcat(PathName, FileName);
        selected_genes=importdata(filename);
        
        Results=CALISTA_path_main(DATA,INPUTS,Results,'selected_genes',selected_genes);
        assignin ('base','Results',Results)
    end
end


% --- Executes during object creation, after setting all properties.
function graphtextwarning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to graphtextwarning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'ForegroundColor','blue');


% --- Executes during object creation, after setting all properties.
function edgestag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgestag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'ForegroundColor','blue');


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
Results_GUI=evalin('base','Results_GUI');

if get(hObject,'Value')==1
    Results_GUI.TRANSITION.jCBList.selectAll;
    %     for k=1:size(h.Results.cluster_distance,1)
    %         jCBList.addCheckBoxListSelectedIndex(k-1);
    %     end
else
    Results_GUI.TRANSITION.jCBList.selectNone;
end
hand = handle(Results_GUI.TRANSITION.jCBList, 'CallbackProperties');
mouseFcn(hand);
assignin ('base','Results_GUI',Results_GUI)

function mouseFcn(hand, ~)
Results_GUI=evalin('base','Results_GUI');
handles_pushbutton1=evalin ('base','handles_pushbutton1');
handles_graphtextwarning=evalin ('base','handles_graphtextwarning');
handles_axes1=evalin('base','handles_axes1');

axes(handles_axes1)
%    fprintf('Item %g selected\n', h.getSelectedIndex+1);
% % Update some items' state programmatically: check/uncheck also selected items
% selected_item=hand.getSelectedIndex;
% 
% if ~isempty(selected_item)
%     
%     temp_checked_items=hand.getCheckBoxListSelectedIndices;
%     if ismember(selected_item,temp_checked_items) %uncheck the selected item if it's checked
%         Results_GUI.jCBList.removeElementAt(selected_item);
%     else
%         Results_GUI.jCBList.addCheckBoxListSelectedIndex(selected_item);
%     end
%     hand = handle(Results_GUI.jCBList, 'CallbackProperties');
% end


fprintf('Item(s) ');
checked_boxes=hand.getCheckBoxListSelectedIndices+1;
fprintf('%g ', checked_boxes);
fprintf('checked\n');
% arrayfun(@cla,findall(0,'type','axes'))
% Update graph
cla(handles_axes1);

if Results_GUI.TRANSITION.transition_new==0
    for k=1:Results_GUI.expected_clusters
        scatter3(Results_GUI.cell_cluster_progression(Results_GUI.final_groups==k)',Results_GUI.score3(Results_GUI.final_groups==k,1),Results_GUI.score3(Results_GUI.final_groups==k,2),30,Results_GUI.colorMARK_calista(k,:),'fill')
        hold on
        text(Results_GUI.TRANSITION.x_center(k),Results_GUI.TRANSITION.y_center(k),Results_GUI.TRANSITION.z_center(k),num2str(k),'FontSize',50);
    end
    
    legend(Results_GUI.legendInfo_calista,'Location', 'northeast')
    xlabel('Cluster pseudotime')
    ylabel('COMP1')
    zlabel('COMP2')
    grid on
    h=[];
    [h]=AddGraph_GUI(Results_GUI.TRANSITION.cluster_distance,Results_GUI.colorMARK_calista,Results_GUI.TRANSITION.x_center,Results_GUI.TRANSITION.y_center,Results_GUI.TRANSITION.z_center,checked_boxes);
else
    h=[];
    [h]=AddGraph_GUI_new(Results_GUI,checked_boxes);

   
   
end
if isempty(h)
    %     Results.TRANSITION.final_graph=h;
    %     h_temp=h;
    %     assignin ('base','h_temp',h_temp)
    set(handles_pushbutton1,'Enable','off')
    set(handles_graphtextwarning,'String','The graph is NOT connected. Please add edges','ForegroundColor','r')
    fprintf('The graph does not exist\n')
else
    NumberOfConnectedNodes=length(dfsearch(h,1));
    
    if NumberOfConnectedNodes==size(h.Nodes,1)
        enableString = get(handles_pushbutton1, 'Enable');
        isEnabled = strcmp(lower(enableString), 'off');
        
        if isEnabled
            set(handles_pushbutton1,'Enable','on')
            set(handles_graphtextwarning,'String','The graph is connected.','ForegroundColor','b')
                 
        end
    else
        enableString = get(handles_pushbutton1, 'Enable');
        isEnabled = strcmp(lower(enableString), 'on');
        if isEnabled
            set(handles_pushbutton1,'Enable','off')
            set(handles_graphtextwarning,'String','The graph is NOT connected. Please add edges','ForegroundColor','r')
            
        end
    end
    
    
end
Results_GUI.TRANSITION.final_graph=h;       
assignin ('base','Results_GUI',Results_GUI)
assignin ('base','handles_pushbutton1',handles_pushbutton1)
assignin ('base','handles_graphtextwarning',handles_graphtextwarning)


function [Results_GUI]=cluster_distance_GUI_update(Results_GUI,handles_axes1)

jCBList=Results_GUI.TRANSITION.jCBList;
% set up callback
hand = handle(jCBList, 'CallbackProperties');
% set(hand, 'MouseClickedCallback', @mouseFcn);
set(hand, 'MousePressedCallback', @mouseFcn);
assignin ('base','Results_GUI',Results_GUI)

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

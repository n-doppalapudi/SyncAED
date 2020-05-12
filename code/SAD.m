function varargout = SAD(varargin)
% SAD MATLAB code for SAD.fig
%      SAD, by itself, creates a new SAD or raises the existing
%      singleton*.
%
%      H = SAD returns the handle to a new SAD or the handle to
%      the existing singleton*.
%
%      SAD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAD.M with the given input arguments.
%
%      SAD('Property','Value',...) creates a new SAD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SAD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SAD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SAD

% Last Modified by GUIDE v2.5 28-Apr-2017 15:38:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SAD_OpeningFcn, ...
                   'gui_OutputFcn',  @SAD_OutputFcn, ...
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


% --- Executes just before SAD is made visible.
function SAD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SAD (see VARARGIN)

% Choose default command line output for SAD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SAD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SAD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function inputfilename_Callback(hObject, eventdata, handles)
% hObject    handle to inputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputfilename as text
%        str2double(get(hObject,'String')) returns contents of inputfilename as a double


% --- Executes during object creation, after setting all properties.
function inputfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inputfilebrowse.
function inputfilebrowse_Callback(hObject, eventdata, handles)
% hObject    handle to inputfilebrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try
    watchon;
% init status
set(handles.StatusText,'String','');

set(handles.StatusText,'String','PMU data is being loaded');


%set(findobj('style','popupmenu'), 'String', {'Date and Time'});
set(handles.inputfilename,'String','filename');
set(handles.outputfilename,'String','Output filename');
%set(handles.plotbrowse,'String','Select Group');

set(handles.totime,'Value',1);
set(handles.totime,'String','Date & Time');

set(handles.fromtime,'Value',1);
set(handles.fromtime,'String','Date & Time');
%set(handles.windowsize,'String','Choose Window size');
%set(findobj('style','popupmenu'), 'String', {'Date and Time'});
%set(findobj('style','listbox'), 'String', {''}, 'Min', 0, 'Max', 1, 'Value', 1, 'ListBoxTop', 1);
set(handles.group,'String','');

%cla(handles.zipplot,'reset');%clf(handles.my_axis_handle);%

global window_size groupno;

global PathEx;
if isempty(PathEx)==1 || PathEx(1) == 0
    PathEx = pwd;
end
[fileName,PathEx] = uigetfile({'*.xlsx;*.csv'},'Select the PMU DATA File',PathEx);
ExPath = [PathEx fileName];
%PathEx(end)=[];

if ExPath(1)==0;%strcmp(ExPath,'null')==1
    
  watchoff;
  msgbox({'File Not Selected'});
  set(handles.StatusText,'String','PMU data was not loaded');
  
  return;
   
else
    [path,name, ext] = fileparts(ExPath);
    guidata(hObject,handles);
    set(handles.inputfilename,'String',ExPath); %% Your editable text box's tag should be "ExLoc" or you should change...
  
end

setPopupmenuString(ExPath, ext, handles);
guidata(hObject,handles);

% windowchoices={60;90;120;150;180;210;240;270;300;};
% set (handles.windowsize,'Value',1)
% set (handles.windowsize,'String',windowchoices)

window_size=60;%% default value if the user does noyt select anything
groupno=1;
guidata(hObject,handles)

set(handles.runanomaly, 'Enable', 'On');

% update status
set(handles.StatusText,'String',sprintf('PMU data of %s is loaded.',fileName));


%function setPopupmenuString(hObject,eventdata,handles)
function setPopupmenuString(hObject,ext,handles)
%fileName=FileEx;
global numbers groupname timeStampString timeStamp;
watchon;
if strcmp(ext,'.xlsx')==1    % excel file
    [numbers,txt] = xlsread(hObject);

    if isempty(txt)
         msgbox({'Input File is not in correct format'});
    return
    end
    
    colNames=txt(1,:)';
   
    timeStampString =txt(2:end,1);% datestr(timeStamp,'yyyy-mm-dd HH:MM:SS.FFF');
    
     if isempty(timeStampString)
         msgbox({'Time Stamp is not in correct format'});
         
         watchoff;
         
     return
     end
    
    timeStamp = datenum(timeStampString);
    
else % csv file   
    
  
    dataTable = readtable(hObject);%, 'HeaderLines',1); % skip 1st line and read data into a table
   
    colNames = dataTable.Properties.VariableNames; % get column names
    %colNames(3) = []; % ingnore status column

    numbers = table2array(dataTable(:,2:end));
    
  
   timeStampString   = table2array(dataTable(:,1));
   
     if isempty(timeStampString)
         msgbox({'Time Stamp is not in correct format'});
         watchoff;
         
     return
     end

   timeStamp=datenum(timeStampString);% = datestr(timeStamp,'yyyy-mm-dd HH:MM:SS.FFF'); 
end

global from_time to_time;

set(handles.fromtime,'Value',1);
%set(handles.fromtime,'string',numbers(:,1));
set(handles.fromtime,'string',timeStampString);

from_time=1;%txt(2:end,1);%1; %% default value if the user does not select anything

set(handles.totime,'Value',length(numbers));
%set(handles.totime,'string',numbers(:,1));
set(handles.totime,'string',timeStampString);

to_time=length(numbers);%txt(end,1); %% default value if the user does noyt select anything
%set(handles.group,'value',numbers)
no_of_groups=(length(colNames)-1)/4;
k=2;
for i=1:no_of_groups
    group_names(i)=strtok(colNames(k),'VM');
    k=k+4;
end

if no_of_groups<1
    msgbox({'Input File is not in correct format'});
    watchoff;
    return
else
    
    
windowchoices={'Replace Bad Data', 'Flag Bad Data'};
set (handles.cleandata,'Value',1)
set (handles.cleandata,'String',windowchoices)
    
groupname=group_names(1);    
set(handles.group,'string',group_names);
end
watchoff;



% --- Executes on selection change in group.
function group_Callback(hObject, eventdata, handles)
% hObject    handle to group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns group contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group
global groupno groupname;
% plotbrowse(get (hObject,'String'));
% guidata(hObject,handles);

groupno=get (hObject,'Value');
guidata(hObject,handles);

contents = get(handles.group,'String'); 
groupname=contents{get(handles.group,'Value')};
guidata(hObject,handles);
%str=handles.plotbrowse;
%set(handles.plotbrowse,'String');
plotname=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function group_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runanomaly.
function runanomaly_Callback(hObject, eventdata, handles)
% hObject    handle to runanomaly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
watchon
% init status
set(handles.StatusText,'String','');
cla(handles.plottype,'reset')
global from_time to_time value cleaning_type numbers groupno untitled outputfile Untitled baddata_detected;

%handles = guidata(hObject);
%set(handles.plotbrowse,'String',groupname); 
guidata(hObject, handles);
%handles.computezip=get(handles.fromtime);
Untitled=numbers((from_time:to_time),(1+(groupno-1)*4:(groupno)*4));
 Gtruth=0;%(V((inserted_baddata),1));
 time=(from_time:to_time)';
 %Untitled=[ Untitled];
untitled= [time Untitled(:,1:2)];
[baddata_index]=MLEprony(untitled,Gtruth);
baddata_detected=zeros(length(time),1);
baddata_detected(baddata_index,1)=1;
%baddata_detected=timeStampString(baddata_index);

% set(handles.outputfilebrowse, 'Enable', 'On');
% guidata(hObject, handles);
% 
% set(handles.saveplot, 'Enable', 'On');
% guidata(hObject, handles);
% 
updatePlot(baddata_detected)




if cleaning_type==1
       for k=1:length(baddata_index)
       indx1=baddata_index(k)+1;
       indx2=baddata_index(k)-1;
       Untitled(baddata_index(k),1)=(Untitled(indx1,1)+Untitled(indx2,1))/2;
       Untitled(baddata_index(k),2)=(Untitled(indx1,2)+Untitled(indx2,2))/2;
       Untitled(baddata_index(k),3)=(Untitled(indx1,3)+Untitled(indx2,3))/2;
       Untitled(baddata_index(k),4)=(Untitled(indx1,4)+Untitled(indx2,4))/2;
       end
      
     outputfile=Untitled;

  plottypes={'Anomaly Plot', 'Original Data Plot', 'Clean Data Plot'};
set (handles.plottype,'Value',1)
set (handles.plottype,'String',plottypes)
% value=get (handles.plottype,'Value');
guidata(hObject, handles);
% handles.plotagain = str2double(plottypes{value});
else if cleaning_type==2
        outputfile=[Untitled baddata_detected];
        plottypes={'Anomaly Plot', 'Original Data Plot'};
set (handles.plottype,'Value',1)
set (handles.plottype,'String',plottypes)
guidata(hObject, handles);
% handles.plotagain = str2double(plottypes{value});
    end
end


%time=(from_time:1:to_time);%timeStampString(from_time:1:to_time);%
set (handles.totime,'Value',to_time);


%outputfile=[a' b' c']; % check this
guidata(hObject, handles);
% update status
set(handles.StatusText,'String','Anomaly Detection is complete');


function updatePlot(baddata_detected)
global from_time to_time timeStamp value;
xData = timeStamp(from_time:to_time);
if value==1
SAD(plot(xData,baddata_detected,'bd'))
title('Anomaly Detection Plot')
ylabel('Anomaly Instance')
elseif value==2
    SAD(plot(xData,baddata_detected))
title('Raw Data Plot')
ylabel('Voltage')
elseif value==3
SAD(plot(xData,baddata_detected))
title('Clean Data Plot')
ylabel('Voltage')
end

datetick('x','hh:MM:ss','keepticks');

ax = gca;
set(ax,'XTickLabelRotation',90) ;
set(ax,'XTickMode', 'auto');

xlabel('Time')


dcm_obj = datacursormode(SAD);
set(dcm_obj,'UpdateFcn',@myupdate_fcn)

if max(baddata_detected)~=0
legend('Anomaly Instances')
end
watchoff;





function outputfilename_Callback(hObject, eventdata, handles)
% hObject    handle to outputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputfilename as text
%        str2double(get(hObject,'String')) returns contents of outputfilename as a double


% --- Executes during object creation, after setting all properties.
function outputfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputfilebrowse.
function outputfilebrowse_Callback(hObject, eventdata, handles)
% hObject    handle to outputfilebrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
watchon;
% init status
set(handles.StatusText,'String','Output file is being written in .xlsx format');

global outputfile groupname timeStampString from_time to_time baddata_detected cleaning_type;
startingFolder = userpath;
defaultFileName = fullfile(startingFolder, '.xlsx');
[baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
if baseFileName == 0
	% User clicked the Cancel button.
    watchoff;
    set(handles.StatusText,'String','User did not select an output file name');

	return;
end
fullFileName = fullfile(folder, baseFileName);
%baseFileName = strcat('''', baseFileName, '''');

if cleaning_type==1

headers(1,1)=groupname;
headers(1,2:5)={'V','VA','I','IA'};
[status, message]=xlswrite(fullFileName,headers);
if status==0
        watchoff;
  msgbox({message});
  set(handles.StatusText,'String',message);
  
  return
end
xlswrite(fullFileName,timeStampString(from_time:to_time),1,'A2')
xlswrite(fullFileName,outputfile,1,'B2')
else if cleaning_type==2
headers(1,1)=groupname;
headers(1,2:6)={'V','VA','I','IA','Bad Data Flag'};
[status, message]=xlswrite(fullFileName,headers);
if status==0
        watchoff;
  msgbox({message});
  set(handles.StatusText,'String',message);
  
  return
end
xlswrite(fullFileName,timeStampString(from_time:to_time),1,'A2')
xlswrite(fullFileName,outputfile,1,'B2')
xlswrite(fullFileName,baddata_detected,1,'F2')
    end
end

set(handles.outputfilename,'String',fullFileName);
watchoff;
% update status
set(handles.StatusText,'String','Output result is saved in a .xlsx file.');

% --- Executes on selection change in cleandata.
function cleandata_Callback(hObject, eventdata, handles)
% hObject    handle to cleandata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cleandata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cleandata
global cleaning_type;
%index_selected = get(hObject,'Value');
% list = (get(hObject,'String'));
% handles.cleandata = str2double(list{index_selected}); 
cleaning_type= get(hObject,'Value');%handles.cleandata;

% --- Executes during object creation, after setting all properties.
function cleandata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cleandata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fromtime.
function fromtime_Callback(hObject, eventdata, handles)
% hObject    handle to fromtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fromtime contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fromtime


% --- Executes during object creation, after setting all properties.
function fromtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fromtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in totime.
function totime_Callback(hObject, eventdata, handles)
% hObject    handle to totime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns totime contents as cell array
%        contents{get(hObject,'Value')} returns selected item from totime


% --- Executes during object creation, after setting all properties.
function totime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function StatusText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StatusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in plottype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottype
global value;
value = get(hObject,'Value');
% list = (get(hObject,'String'));
% handles.plottype = str2double(list{index_selected}); 
% value= handles.plottype;

% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotagain.
function plotagain_Callback(hObject, eventdata, handles)
% hObject    handle to plotagain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global untitled Untitled value baddata_detected;
% index_selected = get(hObject,'Value');
% 
% list = (get(hObject,'String'));
if value==1
    updatePlot(baddata_detected)
else if value==2
        updatePlot(untitled(:,2))
    else if value==3
            updatePlot(Untitled(:,1))
        end
    end
end
    


function figHandle = watchon

%WATCHON Sets the current figure pointer to the watch.
%   figHandle = WATCHON will set the current figure's pointer
%   to a watch.
%
%   See also WATCHOFF.

%   Ned Gulley, 6-21-93
%   Copyright 1984-2014 The MathWorks, Inc.

% If there are no windows open, just set figHandle to a flag value.
if isempty(get(0,'Children')),
    figHandle = NaN;
else
    figHandle = gcf;
    set(figHandle,'Pointer','watch');
    drawnow;
end

function watchoff(figHandle)
%WATCHOFF Sets the given figure pointer to the arrow.
%   WATCHOFF(figHandle) will set the figure's pointer
%   to an arrow. If no argument is given, figHandle is taken to
%   be the current figure.
%
%   See also WATCHON.

%   Ned Gulley, 6-21-93
%   Copyright 1984-2014 The MathWorks, Inc.

if nargin<1
    figHandle = gcf;
end

% If watchon is used before a window has been opened, it will set the
% figHandle to the flag [].  In addition it is generally desirable to not
% error if the window has been closed between calls to watchon and
% watchoff.  ishghandle handles both of these cases.

if ishghandle(figHandle)
    set(figHandle,'Pointer','arrow');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on

% --------------------------------------------------------------------
function Datacursor_Callback(hObject, eventdata, handles)
% hObject    handle to Datacursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacursormode on


% --------------------------------------------------------------------
function Resetinterface_Callback(hObject, eventdata, handles)
% hObject    handle to Resetinterface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StatusText,'String','User has reset the interface');
set(handles.inputfilename,'String','filename');
set(handles.outputfilename,'String','Output filename');
%set(handles.plotbrowse,'String','Select Group');
set(handles.fromtime,'Value',1);
set(handles.totime,'Value',1);
set(handles.fromtime,'String','Date & Time');
set(handles.totime,'String','Date & Time');

set(handles.plottype,'Value',1);
set(handles.plottype,'String','Plot Type');

% set(handles.windowsize,'Value',1);
% set(handles.windowsize,'String','Choose Window size');
set(handles.group,'String','');
set(handles.cleandata,'Value',1);
set(handles.cleandata,'String','Clean Data');
%set(findobj('style','popupmenu'), 'String', {'Date and Time'});
%set('style','listbox', 'String', {''}, 'Min', 0, 'Max', 1, 'Value', 1, 'ListBoxTop', 1);

%set(findobj('style','popupmenu'), 'String', {'Date and Time'});
% set(findobj(handles.inputfilename),'String',{'filename'});
% set(findobj(handles.outputfilename),'String',{'Output filename'});
% set(findobj(handles.plotbrowse),'String',{'Select Group'});
% set(findobj(handles.totime),'String',{'Date & Time'});
% set(findobj(handles.fromtime),'String',{'Date & Time'});
% set(findobj(handles.windowsize),'String',{'Choose Window size'});
% %set(findobj('style','popupmenu'), 'String', {'Date and Time'});
% set(findobj('style','listbox'), 'String', {''}, 'Min', 0, 'Max', 1, 'Value', 1, 'ListBoxTop', 1);
cla(handles.baddataplot,'reset');%clf(handles.my_axis_handle);%
set(handles.runanomaly, 'Enable', 'Off');
set(handles.outputfilebrowse, 'Enable', 'Off');
% set(handles.saveplot, 'Enable', 'Off');
watchoff;

% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;

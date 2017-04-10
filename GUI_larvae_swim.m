function varargout = GUI_larvae_swim(varargin)
% GUI_LARVAE_SWIM MATLAB code for GUI_larvae_swim.fig
%      GUI_LARVAE_SWIM, by itself, creates a new GUI_LARVAE_SWIM or raises the existing
%      singleton*.
%
%      H = GUI_LARVAE_SWIM returns the handle to a new GUI_LARVAE_SWIM or the handle to
%      the existing singleton*.
%
%      GUI_LARVAE_SWIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_LARVAE_SWIM.M with the given input arguments.
%
%      GUI_LARVAE_SWIM('Property','Value',...) creates a new GUI_LARVAE_SWIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_larvae_swim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_larvae_swim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_larvae_swim

% Last Modified by GUIDE v2.5 25-Mar-2016 11:23:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_larvae_swim_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_larvae_swim_OutputFcn, ...
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
end


% --- Executes just before GUI_larvae_swim is made visible.
function GUI_larvae_swim_OpeningFcn(hObject, eventdata, handles, varargin)
clc;cla(handles.picture,'reset');
load iconmov2ima.mat
set(handles.loadimgbutton, 'cdata',Icon);
%%Setappdata
setappdata(0,'hgui',gcf);
put ('toggler',0);
%cla(handles.picture,'reset')
handles.status=0;
handles.output = hObject;
set(handles.picture,'Visible','off');
guidata(hObject, handles);% Update handles structure
end

function varargout = GUI_larvae_swim_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
cla
end

%%
function switchui (who)
handles=guihandles(getappdata(0,'hgui'));
turnoff=findobj('-regexp','Tag','multip');
set(turnoff, 'visible', 'off');
turnon=findobj('-regexp','Tag',who);
set(turnon, 'visible', 'on');
drawnow;
end

function loadimgsbutton_Callback(hObject, eventdata, handles)
[filename, filepath]= uigetfile( ...
    {  '*.jpg','Jpg-files (*.jpg)';'*.jpeg','Jpeg files (*.jpeg)'; ...
    '*.*',  'All Files (*.*)'},'Select the files','MultiSelect', 'on');
filename=filename';
if filepath==0 %if the user pressed cancelled, then we exit this callback
    return
else
    if length(filename)>=1
        m=msgbox('Please wait, loading files...','Loading');
        set(handles.folder_path,'string',filepath);
        
        %% Read images
        % Preallocate the array
        % Calculate number of frames
        nFrames=length(filename);
        I=imread(fullfile(filepath,filename{1}));
        OrigImages=zeros([size(I) nFrames],class(I));%prealocate memory
        %OrigImages(:,:,1)=rgb2gray(I);
        OrigImages(:,:,:,1)=(I);
        % Create image sequence array
        for p=2:nFrames
            %OrigImages(:,:,p)=rgb2gray(imread(fullfile(filepath,filename{p})));
            OrigImages(:,:,:,p)=(imread(fullfile(filepath,filename{p})));
        end
        %save('OrigImages')
        close(m)
        
        % Display picture
        set(handles.picture,'Visible','on');
        handles.nFrames=length(filename);
        handles.OrigImages=OrigImages;
        handles.filepath=filepath;
        handles.filename=filename;
        sliderrange(handles)
        %% Display just first picture
        currentimage=imread([filepath filename{1}]);
        image(currentimage, 'parent',gca, 'cdatamapping', 'scaled');
        colormap('gray');
        vectorcolor='g';
        axis image;
        set(gca,'ytick',[])
        set(gca,'xtick',[])
        set (handles.filenameshow, 'string', ['Frame (' int2str(floor(get(handles.fileselector, 'value'))) '/' int2str(size(filepath,1)/2) '):' sprintf('\n') filename{1}]);
        set (handles.filenameshow, 'tooltipstring', [filepath filename{1}]);
        if strmatch(get(handles.multip01, 'visible'), 'on');
            set(handles.imsize, 'string', ['Image size: ' int2str(size(currentimage,2)) '*' int2str(size(currentimage,1)) 'px' ])
        end
        
    end
    
end %if user press cancel
set (handles.filenamebox, 'string', filename);%displays picture names in popuplist
handles.status=1;%load images
guidata(hObject, handles);% Update handles structure
end

function folder_path_Callback(~, eventdata, handles)
end

function folder_path_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Crop_button_Callback(hObject, eventdata, handles)
%%
h=figure;
imshow(imread(fullfile(handles.filepath,handles.filename{1})),[],'InitialMagnification','fit');
set(h,'color',[1 1 1],'position',[ 1 1 1920 1080]);
title('Click-drag-release to select region of interest')
% user select region
done=false;
while ~done
    text(-1740.9,-320,'Please take advantage to double check that the camera is alineated with the column plane','FontSize',14)
    [subIm,roi]=imcrop;
    %wait for user to select region of interest (click-drag rubberband box)
    %determine corners of selected region
    x1=roi(1);  x2=roi(1)+roi(3);  xx=[x1 x1 x2 x2 x1];
    y1=roi(2);  y2=roi(2)+roi(4);  yy=[y1 y2 y2 y1 y1];
    %draw rectangle to indicate selected region
    hLine=line(xx,yy);
    title('Click inside to keep, outside to select again')
    pt=ginput(1);
    %wait for user to keep or discard selected region
    % check: is point inside region?
    if pt(1)>=roi(1) && pt(1)<=roi(1)+roi(3) && pt(2)>=roi(2) && pt(2)<=roi(2)+roi(4)
        done=true;  %exit loop
    else
        delete(hLine)   %remove previous selection
        title('Click-drag-release to select region of interest')
        %loop around to select another region
    end
end
hold on
plot(x1,y2,'ro')
text(x1*2,y2,' Reference point 0,0')
title('Done')
Xpixels=x1-roi(1);Ypixels=y2-roi(2);
%pause(1)
close(h)
%% =======================================================================
%crop all frames; display cropped footage
for i=1:handles.nFrames
    Crop_OrigImages(:,:,:,i)=imcrop(handles.OrigImages(:,:,:,i),roi);
    Images(:,:,i)=rgb2gray(Crop_OrigImages(:,:,:,i));
end
clear done hLine pt x1 x2 xx y1 y2 yy %subIm
%We Will use this at the end
%CropImg=Images(:,:,1)+Images(:,:,2);
%% Storage and preallocate data and structures
handles.Crop_OrigImages=Crop_OrigImages;
handles.Images=Images;
handles.bw2=zeros(size(Images));
handles.subIm=subIm;
centroid.x=[];centroid.y=[];centroid = repmat(centroid,handles.nFrames,1);
handles.centroid=centroid;
handles.roi=roi;
%%========================================================================
%% display croped Images
%imshow(Images(:,:,str2double(get(handles.Frame_No,'String'))),[],'InitialMagnification','fit')
sliderdisp(hObject,eventdata,handles)
%save('Images')
%%========================================================================
handles.status=2;%crop images images
guidata(hObject, handles);% Update handles structure
hold off
end

function Button_particle_segmentation_Callback(hObject, eventdata, handles)
Message= msgbox('Please wait','Processing image','help');
Images=handles.Images;
%handles.bw2=[];
%%
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
a=Images(:,:,selected);handles.a=a;
%gfr=imcomplement(a);
%imshow(gfr,[],'InitialMagnification','fit')
background = imopen(imcomplement(a),strel('disk',str2double(get(handles.radius,'String'))));
bw = im2bw(imadjust(imcomplement(a)- background),str2double(get(handles.threshold,'String')));
[~,L]=bwboundaries(bw(:,:,1),'noholes');
stats=regionprops(L,'area','PixelIdxList');%Measure properties of image regions (e.g Area
id=[stats.Area]>str2double(get(handles.P,'String'));
stats_id=stats(id);
idxToKeep=vertcat(stats_id(:).PixelIdxList);
bw2= false(size(bw));
bw2(idxToKeep) = true;
%imshow(bw2,[],'InitialMagnification','fit')
bw2=imclose(bw2,strel('disk',str2double(get(handles.gaps,'String'))));
handles.bw2(:,:,selected)=bw2;
set(handles.Identify_centroid_button,'Visible','on');
guidata(hObject, handles);% Update handles structure
delete(Message)
sliderdisp(hObject,eventdata,handles)
handles.status=3;%Preview image in binary
guidata(hObject, handles);% Update handles structure
end
% --------------------------------------------------------------------

function Convert_all_Frames2binary_menu_Callback(hObject,~, handles)

%Message= msgbox('Please wait','Processing images','help');
%subIm=handles.subIm;
Images=handles.Images;
%Images3=false([size(subIm,1) size(subIm,2) handles.nFrames]);
%handles.bw2=[];
%handles.a=[];
%% Apply to all
h = waitbar(0,'Converting images to binary form...','Name','Larvae detection...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
if getappdata(h,'canceling')
    delete(h);
    Exit=1;
    return;
end
%waitstep = floor((handles.nFrames)/100);%if there are more than 100 pictures
%set(handles.text_Processing_frame,'visible','on');
%title(['Time=',num2str(round((time(round(t/dt_multip05+1))/3600)*10)/10),' hours','   ','Egg diameter=',num2str(round(D(round(t/dt_multip05+1))*10)/10),' mm'],'FontSize',fontsize)
%% ========================================================================
fill=0.4/(handles.nFrames+0.4); %fills waitbar with time spent until this point
% Check for Cancel button press
if getappdata(h,'canceling')
    delete(h);
    Exit=1;
    return;
end
% Report current estimate in the waitbar's message field
waitbar(fill,h,['Please wait....' sprintf('%12.0f',fill*100) '%']);
%% ======================================
for i=1:handles.nFrames
    %%
    %%
    %if sum(sum(handles.bw2(:,:,i)))==0
    %set(handles.text_Processing_frame,'string',['Processing frame ',num2str(i),' of ',num2str(handles.nFrames)])
    %a=Images(:,:,i);
    gfr=imcomplement(Images(:,:,i));
    %imshow(gfr,[],'InitialMagnification','fit');
    background = imopen(gfr,strel('disk',str2double(get(handles.radius,'String'))));
    bw = im2bw(imadjust(gfr - background),str2double(get(handles.threshold,'String')));
    [~,L]=bwboundaries(bw(:,:,1),'noholes');
    stats=regionprops(L,'area','PixelIdxList');%Measure properties of image regions (e.g Area
    id=[stats.Area]>str2double(get(handles.P,'String'));
    stats_id=stats(id);
    idxToKeep=vertcat(stats_id(:).PixelIdxList);
    bw2= false(size(bw));
    bw2(idxToKeep) = true;
    handles.bw2(:,:,i)=imclose(bw2,strel('disk',str2double(get(handles.gaps,'String'))));
    %end
    %if ~mod(i, waitstep) || i==handles.nFrames %if there are more than 100 pictures update the waitbar periodically
    fill=(i)/(handles.nFrames+0.4);%0.4 is the approx initial time before the loop
    % Check for Cancel button press
    if getappdata(h,'canceling')
        delete(h);
        Exit=1;
        return;
    end
    % Report current estimate in the waitbar's message field
    waitbar(fill,h,['Please wait....' sprintf('%12.0f',fill*100) '%']);
    %end%if mod
end%for

%% DELETE the waitbar;
delete(h)
%bw2=Images3;

%handles.bw2=bw2;
%handles.a=a;
%%
%display('Done')
set(handles.Plot_all_button,'Visible','on');
% set(handles.Centroid_all,'Visible','on');
% set(handles.text_timeinterval_multip05,'Visible','on');
% set(handles.dt_multip05,'Visible','on');
%delete(Message)
handles.status=4;%Convert all to image in binary
guidata(hObject, handles);% Update handles structure
guidata(hObject, handles);% Update handles structure
imshow(bw2(:,:,end),[],'InitialMagnification','fit');
Message= msgbox('Done','Done','help');
end
function Identify_centroid_button_Callback(hObject, eventdata, handles)
%% Load data from handles
bw2=handles.bw2;%a=handles.a;
toggler=retr('toggler');

%% Which frame you want to display??
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);

%Check if the image was converted to binary
if sum(sum(bw2(:,:,selected)))==0
    errordlg('Please convert the image to binary form','Error');
    return
end

[~,L]=bwboundaries(bw2(:,:,selected),'noholes');
stats2=regionprops(L,'centroid');%Measue properties of image regions (e.g Area
x2=[stats2.Centroid]';x2=x2(1:2:end);%get just the x coordinate
y2=[stats2.Centroid]';y2=y2(2:2:end);%get just the y coordinate
%% Figure display
Orig_croped=imcrop(imread(fullfile(handles.filepath,handles.filename{selected})),handles.roi);
%h=figure;set(h,'color',[1 1 1],'position',[ 1 1 1920 1080]);
%subaxis(1,2,1,'mt',0.01,'mb',0.04);
imshow(Orig_croped,[],'InitialMagnification','fit');
hold on;plot(x2,y2,'ro');
%subaxis(1,2,2,'mt',0.01,'mb',0.04);
imshow(bw2(:,:,selected),[],'InitialMagnification','fit');hold on;plot(x2,y2,'ro')
%% Storage data
handles.centroid(selected, 1).x=x2;
handles.centroid(selected, 1).y=y2;
guidata(hObject, handles);% Update handles structure
hold off
end
function P_Callback(hObject, eventdata, handles)
end
function P_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function gaps_Callback(hObject, eventdata, handles)
end
function gaps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function Centroid_all_Callback(hObject, eventdata, handles)
bw2=handles.bw2;Images=handles.Images;
% centroid.x=[];centroid.y=[];
% centroid = repmat(centroid,handles.nFrames,1);
%i=str2double(get(handles.Frame_No,'String'));
%figure('color','w');imshow(bw2(:,:,1),[],'InitialMagnification','fit','Border','tight')
for i=1:handles.nFrames
    if isempty(handles.centroid(i, 1).x)%If we have not identified the centroid do this...
        [~,L]=bwboundaries(bw2(:,:,i),'noholes');
        stats2=regionprops(L,'all');%Measure properties of image regions (e.g Area
        x2=[stats2.Centroid]';x2=x2(1:2:end);%get x coordinates
        y2=[stats2.Centroid]';y2=y2(2:2:end);%get y coordinates
        handles.centroid(i, 1).x=x2; handles.centroid(i, 1).y=y2;%storage data in centroid structure
    end
    if i==handles.nFrames
        set(handles.fileselector, 'value',handles.nFrames)
        plot_Image(handles,bw2(:,:,i))
        hold on;plot(handles.centroid(i, 1).x,handles.centroid(i, 1).y,'ro')
    end
end


%handles.centroid=centroid;
Message= msgbox('Done','Done','help');
handles.status=5;%Identify cetroid in all
guidata(hObject, handles);% Update handles structure
hold off
end

% --- Executes on button press in Track_larvae_button.
function Track_larvae_button_Callback(hObject, eventdata, handles)
nFrames=handles.nFrames;

if strcmp(get(handles.multip06, 'visible'), 'on')
    %% ==================================================================================
    %% Storage data
    % Create points cell
    %
    cla reset
    max_linking_distance =str2double(get(handles.max_linking_distance,'string')) ;%1000;%343.7;
    max_gap_closing = str2double(get(handles.max_gap_closing,'string'));%Inf
    debug = true;
    points = cell(nFrames, 1);times=points;%allocate structures
    %%=== TIMES ===========================================================
    Frames_time=Detect_CameraDateTime(handles);
    %Put data in points and times structures
    for i=1:nFrames
        points{i} = [(handles.centroid(i, 1).x) (handles.centroid(i, 1).y)];
        times{i}=Frames_time(i)*ones(length(handles.centroid(i, 1).x),1);
    end
    %
    [ tracks adjacency_tracks ] = simpletracker(points,...
        'MaxLinkingDistance', max_linking_distance, ...
        'MaxGapClosing', max_gap_closing, ...
        'Debug', debug,'Method','Hungarian');%'NearestNeighbor'
end %If we are in the module to identify paths

%% Acess data if not in tracking larvae path module
if strcmp(get(handles.multip07, 'visible'), 'on')
    points=handles.points;
    tracks=handles.tracks;
    adjacency_tracks=handles.adjacency_tracks;
    times=handles.times;
end
%% Plot centroid of larvae location
CM = hsv(nFrames); %colormap

for j=1:nFrames
    plot(points{j, 1}(:,1),points{j, 1}(:,end),'MarkerFaceColor',CM(j,:),'MarkerEdgeColor','k','Marker','o','LineStyle','none');
    hold on
    set(gca,'YDir','reverse')
    legendInfo{j} = [handles.filename{j,1}];
end
h=legend(legendInfo);
set(h,'interpreter','none')
n_tracks = numel(tracks);
colors =lines(n_tracks);
all_points = vertcat(points{:});
all_times=vertcat(times{:});
larvae_tracks = cell(n_tracks, 1);
for i_track = 1 : n_tracks  %n_tracks=number of paths
    
    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    time_consecutive_larvae{i_track}=all_times(track, :);
    larvae_tracks{i_track}=[track_points(:,1)  track_points(:, 2)];
    plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :),'linewidth',2)
    %quiver(track_points(:,1), track_points(:,
    %2),6000*ones(size(track_points,1),1),6000*ones(size(track_points,1),1))%for
    %later
    hold all
end

%Look for lines and send them to the back.
lineObj=findobj(gca,'type','line');
%uistack(lineObj, 'down',1)
if strcmp(get(handles.multip06, 'visible'), 'on')
    xlabel('X, in pixels')
    ylabel('Y, in pixels')
end
hold off
handles.larvae_tracks=larvae_tracks;
handles.time_consecutive_larvae=time_consecutive_larvae;
handles.points=points;
handles.track_points=track_points;
handles.tracks=tracks;
handles.adjacency_tracks=adjacency_tracks;
handles.times=times;
handles.status=6;%Track larvae
guidata(hObject, handles);% Update handles structure
%% =======================================================================================
    function [seconds]=Detect_CameraDateTime(handles)
        filepath=handles.filepath;
        filename=handles.filename;
        nFrames=size(filename,1);
        %% Get camera info
        for i=1:nFrames
            info(i)=imfinfo(fullfile(filepath,filename{i}));%imfinfo->gets the file metadata.  fullfile-->generates the full path for every figure
        end
        DigCamera=arrayfun(@(cam) cam.DigitalCamera,info); date={DigCamera.DateTimeOriginal}';
        Milliseconds=str2double({DigCamera.SubsecTimeOriginal}')/100;
        datevect=datevec(date,'yyyy:mm:dd HH:MM:SS');
        datevect(:,end)=datevect(:,end)+Milliseconds;
        seconds=(datevect(:,4)*3600)+(datevect(:,5)*60)+datevect(:,end);
    end%end function Detect_CameraDateTime
end


% --- Executes on button press in Manually_delete_path.
function Manually_delete_path_Callback(hObject, eventdata, handles)
set(handles.Text_messages,'string','Please click on the last larvae centroid of the path you want to delete.') 
[X,Y]=ginput(1);
minDist      = Inf;
minHndlIdx   = 0;

%% Acess data
points=handles.points;
tracks=handles.tracks;
adjacency_tracks=handles.adjacency_tracks;
all_points = vertcat(points{:});

% identify larvae location in points structure

%Use the X and Y data for each larvae at each frame (all points
%structure)and determine which one is closest to the mouse down event
xData =all_points(:,1);
yData =all_points(:,2);
dist  =((xData-X).^2+(yData-Y).^2);
minHndlIdx=find(dist==min(dist),1);
minDist = min(dist);

% In which path is the larvae?
n_tracks = numel(tracks);
larvae_tracks = cell(n_tracks, 1);
for i_track = 1 : n_tracks  %n_tracks=number of paths
    if ~isempty(find(adjacency_tracks{i_track}==minHndlIdx))
        Path_to_modify=i_track;
        tracks{i_track}=[];
        adjacency_tracks{i_track}=[];
    end
end
adjacency_tracks(cellfun(@(adjacency_tracks) isempty(adjacency_tracks),adjacency_tracks))=[];
tracks(cellfun(@(tracks) isempty(tracks),tracks))=[];
% handles.larvae_tracks=larvae_tracks;
% handles.time_consecutive_larvae=time_consecutive_larvae;
% handles.points=points;
% handles.track_points=track_points;
handles.tracks=tracks;
handles.adjacency_tracks=adjacency_tracks;
% handles.times=times;
guidata(hObject, handles);% Update handles structure
Track_larvae_button_Callback(hObject, eventdata, handles)
set(handles.Text_messages,'string',' ') 
end


% --- Executes on button press in add_larvae_path.
function add_larvae_path_Callback(hObject, eventdata, handles)
set(handles.Text_messages,'string','Please click on the centroid of each larvae you want to add to the path and right click in the centroid of the last larvae of the path.') 
x= [];
y = [];
points=handles.points;
all_points = vertcat(points{:});
xData =all_points(:,1);
yData =all_points(:,2);%lavae location
frames=zeros(size(xData));
% counter=1;
% for i=1:length(points)  
%     frames(counter:counter+size(points{i,1},1)-1)=repmat(i,1,size(points{i,1},1))
%     counter=counter+size(points{i,1},1);
% end
% n_slices = numel(points);
done=false;
while ~done
    [xi,yi,button] = ginput(1);
    x = [x;xi];
    y = [y;yi];
    if button==3
        done=true;  %exit
        break
    end 
end
Idtracks=zeros(1,length(x));
for i=1:length(x)   
dist  =((xData-x(i)).^2+(yData-y(i)).^2);
Idtracks(i)=find(dist==min(dist),1);
end
%% Create a path
adjacency_tracks=handles.adjacency_tracks;
adjacency_tracks{end+1}=Idtracks;
tracks=handles.tracks;

n_cells = cellfun(@(x) size(x, 1), points); %number of points in each picture
 track = NaN(numel(points), 1);
 adjacency_track=adjacency_tracks{end};
        
        for j = 1 : numel(adjacency_track)
            
            cell_index = adjacency_track(j);%point number from adjacency number
            
            % We must determine the frame this index belong to
            tmp = cell_index;
            frame_index = 1;
            while tmp > 0
                tmp = tmp - n_cells(frame_index);
                frame_index = frame_index + 1; %find picture where call_index point number is located based on n-cells array
            end
            frame_index = frame_index - 1; 
            in_frame_cell_index = tmp + n_cells(frame_index); %number of points in corresponding picture
            
            track(frame_index) = in_frame_cell_index;
            
        end
        
        tracks{end+1} = track;
        
handles.tracks{end+1}=tracks;
handles.adjacency_tracks=adjacency_tracks;
guidata(hObject, handles);% Update handles structure
Track_larvae_button_Callback(hObject, eventdata, handles)
set(handles.Text_messages,'string',' ') 
end


% % --- Executes on button press in Select_reference_points_button.
function Select_reference_points_button_Callback(hObject, eventdata, handles)
%n=input('number of points to use in the calibration ');
n_sq_x=str2double(get(handles.n_sq_x,'string'));%number of squares
n_sq_y=str2double(get(handles.n_sq_y,'string'));
n=(n_sq_x+1)*(n_sq_y+1);
dx=str2double(get(handles.dx,'string'));%1cm
dy=str2double(get(handles.dy,'string'));%1cm
%
% Annotation=annotation(h,'textbox',[0.45 0.8 0.28 0],...
%     'String',{'Zoom in and click in consecutive points'},...
%     'FontWeight','bold','FontSize',12,'FitBoxToText','off',...
%     'EdgeColor','none','Color',[1 1 1]);
% title('Click on the figure every 1cm');
% disp('Click on the figure every 1cm');
x= [];
y = [];
hold on;
%pause
for count=1:n,
    [xi,yi] = ginput(1);
    plot(xi,yi,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);
    x = [x;xi];
    y = [y;yi];
    if count==1
        Xpixels=xi;
        Ypixels=yi;
    end%end if
end%end for

%pixPerCm=Diff(y)/1;
Xcm=0;Ycm=0;
pixels=[x-Xpixels y-Ypixels];
clear xi yi count
%pixels=[x y];
if n_sq_x==0%case of a ruller
    cm=[zeros(n,1) (0:dy:n-1)'];
else
    y_l = ((0:n_sq_y)'*ones(1,n_sq_x+1)); x_l = (ones(n_sq_y+1,1)*(0:n_sq_x));
    cm= [x_l(:) y_l(:)];
    if (n_sq_y+1)*(n_sq_x+1)>n
        cm(n+1:(n_sq_y+1)*(n_sq_x+1),:)=[];
    end
end
%Least squares pixels=cm*A --> A=[pixeles/cm]
%pixel/cm ratio=lscov(cm,pixels)  -->cm*pixel/cm ratio=pixels
%A=lscov(cm,pixels);%To get cm-->cmx=pixels(:,1)/A(1,1);
Pix_Cmy=lscov(cm(:,2),pixels(:,2));
Pix_Cmx=lscov(cm(:,1),pixels(:,1));
%figure;
plot(cm(:,1)*Pix_Cmx+Xpixels,cm(:,2)*Pix_Cmy+Ypixels,'o')
set(handles.Pixel2cm_ratio_text,'string',['The pixel by cm ratio is equal to ' num2str(round((Pix_Cmx)*10)/10) ' and '  num2str(round((Pix_Cmy)*10)/10) ' in x and y, respectively.' ] )
handles.Pix_Cmx=Pix_Cmx;
handles.Pix_Cmy=Pix_Cmy;
handles.status=7;%Calibrated
guidata(hObject, handles);% Update handles structure
hold off
end

% --------------------------------------------------------------------
function Calculate_velocities_Callback(hObject, eventdata, handles)
set(handles.tools,'visible','off')
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip09')
if ~isfield(handles,'Pix_Cmx')
    errordlg('Please process the images and calibrate camera','Error');
    return
end
Pix_Cmx=handles.Pix_Cmx;
Pix_Cmy=handles.Pix_Cmy;
larvae_tracks=handles.larvae_tracks;
time_consecutive_larvae=handles.time_consecutive_larvae';
Larvae_velocity_x_and_y=cell(length(larvae_tracks),1);
Distance_x_y=cell(length(larvae_tracks),1);
time_span=cell(length(larvae_tracks),1);
Larvae_velocity_Mag=Larvae_velocity_x_and_y;
mean_vel_x_and_y=[];
mean_Distance_x_y=[];
mean_time_span=[];
larvae_withPath=1; %For larvae that belong to a path
for i=1:length(larvae_tracks)
    if size(larvae_tracks{i},1)>1
        Larvae_velocity_x_and_y{larvae_withPath}=[diff(larvae_tracks{i})./[(diff(time_consecutive_larvae{i}))*Pix_Cmx (diff(time_consecutive_larvae{i}))*Pix_Cmy]];
        Larvae_velocity_Mag{larvae_withPath}=sqrt(Larvae_velocity_x_and_y{larvae_withPath, 1}(:,1).^2)+(Larvae_velocity_x_and_y{larvae_withPath, 1}(:,end).^2);
        Distance_x_y{larvae_withPath}=diff(larvae_tracks{i})/Pix_Cmx;%in cm
        time_span{larvae_withPath}=diff(time_consecutive_larvae{i});%seconds
        if size(Larvae_velocity_x_and_y{larvae_withPath},1)>1
            mean_vel_x_and_y(larvae_withPath,:)= mean(Larvae_velocity_x_and_y{larvae_withPath},1);
            mean_Distance_x_y(larvae_withPath,:)= mean(Distance_x_y{larvae_withPath},1);
            mean_time_span(larvae_withPath,:)= mean(time_span{larvae_withPath},1);
        else
            mean_vel_x_and_y(larvae_withPath,:)= Larvae_velocity_x_and_y{larvae_withPath};
            mean_Distance_x_y(larvae_withPath,:)= Distance_x_y{larvae_withPath};
            mean_time_span(larvae_withPath,:)= time_span{larvae_withPath};
        end
        larvae_withPath=larvae_withPath+1;
    end
end%end for
mean_Vel_mag=sqrt((mean_vel_x_and_y(:,1).^2)+(mean_vel_x_and_y(:,2).^2));
handles.mean_Vel_mag=mean_Vel_mag;
handles.mean_vel_x_and_y=mean_vel_x_and_y;
handles.mean_Distance_x_y=mean_Distance_x_y;
handles.mean_time_span=mean_time_span;
handles.Larvae_velocity_x_and_y=Larvae_velocity_x_and_y;
handles.Distance_x_y=Distance_x_y;
handles.time_span=time_span;

handles.status=8;%Calibrated
guidata(hObject, handles);% Update handles structure
Message= msgbox('Done','Done','help');
pause(0.8)
if exist('Message', 'var')
    delete(Message);
    clear('Message');
end
end%end function

function sliderrange(handles)
%filepath=handles.filepath;
filename=handles.filename;
if size(filename,1)>2
    sliderstepcount=size(filename,1);
    set(handles.fileselector, 'enable', 'on');
    set (handles.fileselector,'value',1, 'min', 1,'max',sliderstepcount,'sliderstep', [1/(sliderstepcount-1) 1/(sliderstepcount-1)*10]);
else
    %sliderstepcount=1;
    set(handles.fileselector, 'enable', 'off');
    set (handles.fileselector,'value',1, 'min', 1,'max',2,'sliderstep', [0.5 0.5]);
end
end
function fileselector_Callback(hObject, eventdata,handles)
filename=handles.filename;
if size(filename,1) > 1
    sliderdisp(hObject,eventdata,handles)
end
end
function sliderdisp(hObject,eventdata,handles)
cla(handles.picture,'reset')
filepath=handles.filepath;
filename=handles.filename;
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
if selected==0 % If return to initial frame while toggle is activated, setup selected to first frame
    selected=1;
end
%% Image to display
if (strcmp(get(handles.multip01, 'visible'), 'on')+0)+(strcmp(get(handles.multip08, 'visible'), 'on')+0)>=1;
    currentimage=imread([filepath filename{selected}]);
elseif (strcmp(get(handles.multip02, 'visible'), 'on')+0)+(strcmp(get(handles.multip07, 'visible'), 'on')+0)>=1;
    Crop_OrigImages=handles.Crop_OrigImages;
    currentimage=Crop_OrigImages(:,:,:,selected);
elseif (strcmp(get(handles.multip03, 'visible'), 'on')+0)+(strcmp(get(handles.multip04, 'visible'), 'on')+0)>=1
    currentimage=handles.bw2(:,:,selected);
elseif strcmp(get(handles.multip05, 'visible'), 'on')
    switch get(handles.popupmenu_BackgroundImage,'value')
        case 1 %'Binary image'
            currentimage=handles.bw2(:,:,selected);
        case 2 %'Cropped original image'
            currentimage=handles.Crop_OrigImages(:,:,:,selected);
    end
    
end
plot_Image(handles,currentimage)
if strcmp(get(handles.multip05, 'visible'), 'on')
    hold on;
    centroid=handles.centroid;
    plot(centroid(selected).x,centroid(selected).y,'ro')
    hold off
end
if strcmp(get(handles.multip07, 'visible'), 'on')
    hold on;
    %% Acess data if not in tracking larvae path module
    points=handles.points;
    tracks=handles.tracks;
    adjacency_tracks=handles.adjacency_tracks;
    times=handles.times;
    Track_larvae_button_Callback(hObject, eventdata, handles)
    hold off
end
%hold on;
%plot(handles.centroid(selected, 1).x,handles.centroid(selected, 1).y,'ro')
%hold off
end
function plot_Image(handles,currentimage)
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
if selected==0 % If return to initial frame while toggle is activated, setup selected to first frame
    selected=1;
end
filepath=handles.filepath;
filename=handles.filename;
%%========================
%image(currentimage, 'parent',gca,'cdatamapping', 'scaled');
imshow(currentimage,[],'InitialMagnification','fit')
% colormap('gray');
% vectorcolor='g';
axis image;
set(gca,'ytick',[])
set(gca,'xtick',[])
if toggler==0
    set (handles.filenameshow, 'string', ['Frame (' int2str(floor(get(handles.fileselector, 'value'))) '/' int2str(size(filepath,1)/2) '):' sprintf('\n') filename{selected}]);
else
    set (handles.filenameshow, 'string', ['Frame (' int2str(selected) '/' int2str(size(filepath,1)/2) '):' sprintf('\n') filename{selected}]);
end
set (handles.filenameshow, 'tooltipstring', [filepath filename{selected}]);
if strmatch(get(handles.multip01, 'visible'), 'on');
    set(handles.imsize, 'string', ['Image size: ' int2str(size(currentimage,2)) '*' int2str(size(currentimage,1)) 'px' ])
end
drawnow;
end
function put(name, what)
hgui=getappdata(0,'hgui');
setappdata(hgui, name, what);
end
function var = retr(name)
hgui=getappdata(0,'hgui');
var=getappdata(hgui, name);
end
function threshold_Callback(hObject, eventdata, handles)
end
function threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function togglepair_Callback(hObject, eventdata, handles)
toggler=get(gco, 'value');
put ('toggler',toggler);
%filepath=handles.filepath;
filename=handles.filename;
%filepath=retr('filepath');
if size(filename,1) > 1
    sliderdisp(hObject,eventdata,handles)
end
end

% --- Executes during object creation, after setting all properties.
function fileselector_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function Preprocessing_menu_Callback(hObject, eventdata, handles)
if isfield(handles,'filename')==0
    errordlg('Please load images','Error');
end
end

function New_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip01')
if strcmp(get(handles.tools, 'visible'), 'off');
    set(handles.tools,'visible','on')
end
end
% --------------------------------------------------------------------
function Select_Region_of_interest_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.fileselector, 'value',1)
currentimage=handles.OrigImages(:,:,:,1);
plot_Image(handles,currentimage)
switchui('multip02')
if strcmp(get(handles.tools, 'visible'), 'off');
    set(handles.tools,'visible','on')
end
end
% --------------------------------------------------------------------
function Convert2binary_settings_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip03')
if strcmp(get(handles.tools, 'visible'), 'off');
    set(handles.tools,'visible','on')
end
end
% --------------------------------------------------------------------
function Binary_all_Frames_menu_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip04')
if strcmp(get(handles.tools, 'visible'), 'off');
    set(handles.tools,'visible','on')
end
Convert_all_Frames2binary_menu_Callback(hObject,eventdata,handles)
end
% --------------------------------------------------------------------
function Detect_Larvae_in_all_frames_menu_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip05')
if strcmp(get(handles.tools, 'visible'), 'off');
    set(handles.tools,'visible','on')
end
Centroid_all_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function Track_Larvae_Paths_menu_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip06')
set(handles.tools,'visible','off')
end

% --------------------------------------------------------------------
function Reject_path_Callback(hObject, eventdata, handles)
switchui('multip07')
set(handles.tools,'visible','on')
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
if selected==0 % If return to initial frame while toggle is activated, setup selected to first frame
    selected=1;
end
%% Image to display
Crop_OrigImages=handles.Crop_OrigImages;
currentimage=Crop_OrigImages(:,:,:,selected);
plot_Image(handles,currentimage)
%% Add on scatter of larvae centroids with larvae paths
hold on
Track_larvae_button_Callback(hObject, eventdata, handles)
hold off
%handles.status=7;%Manually delete larvae
%guidata(hObject, handles);% Update handles structure
end

% --------------------------------------------------------------------
function Calibration_menu_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui('multip08')
set(handles.fileselector, 'value',1)
currentimage=handles.OrigImages(:,:,:,1);
plot_Image(handles,currentimage)
set(handles.tools,'visible','on')
end

% --- Executes on button press in Plot_all_button.
function Plot_all_button_Callback(hObject, eventdata, handles)
bw2=handles.bw2;
%% =======================================================================
for i=1:handles.nFrames
    set(handles.fileselector, 'value',i)
    plot_Image(handles,bw2(:,:,i))
    %title(sprintf('Frame %d',i))
    %drawnow
    pause(str2double(get(handles.dt_multip05,'String')))
end
end


function Plot_all_centroid_button_Callback(hObject, eventdata, handles)
centroid=handles.centroid;
switch get(handles.popupmenu_BackgroundImage,'value')
    case 1 %'Binary image'
        currentimages=handles.bw2;
    case 2 %'Cropped original image'
        currentimages=handles.Crop_OrigImages; %currentimage=handles.Crop_OrigImages(:,:,:,selected);
end

for i=1:handles.nFrames
    cla(handles.picture,'reset')
    %h=figure;set(h,'color',[1 1 1],'position',[ 1 1 1920 1080]);
    %subplot(1,2,1);imshow(Images(:,:,i),[],'InitialMagnification','fit');
    %title(sprintf('Frame %d',i))
    set(handles.fileselector, 'value',i)
    %%
    switch get(handles.popupmenu_BackgroundImage,'value')
        case 1 %'Binary image'
            currentimages=handles.bw2(:,:,i);
        case 2 %'Cropped original image'
            currentimages=handles.Crop_OrigImages(:,:,:,i); %currentimage=handles.Crop_OrigImages(:,:,:,selected);
    end
    %%
    plot_Image(handles,currentimages)
    %title(sprintf('Frame %d',i))
    hold on;
    plot(centroid(i).x,centroid(i).y,'ro')
    hold off
    %     imshow(bw2(:,:,i),[],'InitialMagnification','fit');hold on;plot(x,y,'ro')
    %     title(sprintf('Frame %d',i))
    pause(str2double(get(handles.dt_multip05,'String')))
    drawnow
    
end
set(gcf,'doublebuffer','off')
end

function popupmenu_BackgroundImage_Callback(hObject, eventdata, handles)
centroid=handles.centroid;
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
if selected==0 % If return to initial frame while toggle is activated, setup selected to first frame
    selected=1;
end

%% ========================================================================
switch get(handles.popupmenu_BackgroundImage,'value')
    case 1 %'Binary image'
        currentimage=handles.bw2(:,:,selected);
    case 2 %'Cropped original image'
        currentimage=handles.Crop_OrigImages(:,:,:,selected);
end
plot_Image(handles,currentimage)

%% ========================================================================
hold on
plot(centroid(selected).x,centroid(selected).y,'ro')
hold off

if get(handles.popupmenu_BackgroundImage,'value')==2
    set(handles.multip5a,'visible','on')
end
if get(handles.popupmenu_BackgroundImage,'value')==2
    set(handles.multip5b,'visible','on')
end

end


% --------------------------------------------------------------------
function Load_sesion_Callback(hObject, eventdata, handles)
[Session_FileName,Session_PathName, filterindex] = uigetfile({'*.mat','MATLAB Files (*.mat)'; '*.mat','mat'},'Load LarvaeSwim session');
ed=msgbox('Loading session, please wait','Loading','help');
status=handles.status;
cla(handles.picture,'reset');
if exist('ed', 'var')
    delete(ed);
    clear('ed');
end
if isequal(Session_FileName,0) | isequal(Session_PathName,0)
else
    clear iptPointerManager
    warning off all
    vars=load(fullfile(Session_PathName,Session_FileName));
    names=fieldnames(vars);
    for i=1:size(names,1)
        handles.(names{i})=vars.(names{i});
    end
    guidata(hObject, handles);% Update handles structure
    sliderrange(handles)
end

if status==1
    set(handles.fileselector, 'value',1)
    if strcmp(get(handles.tools, 'visible'), 'off');
        set(handles.tools,'visible','on')
    end
    currentimage=handles.OrigImages(:,:,:,1);
    plot_Image(handles,currentimage)
    %% settings
    set(handles.folder_path,'string',handles.filepath);
    set (handles.filenameshow, 'string', ['Frame (' int2str(1) '/' int2str(size(handles.filepath,1)/2) '):' sprintf('\n') handles.filename{1}]);
    set (handles.filenameshow, 'tooltipstring', [handles.filepath handles.filename{1}]);
    set(handles.imsize, 'string', ['Image size: ' int2str(size(currentimage,2)) '*' int2str(size(currentimage,1)) 'px' ])
    set (handles.filenamebox, 'string', handles.filename);%displays picture names in popuplist
    
elseif status==2 %crop
    set(handles.fileselector, 'value',1)
    currentimage=handles.OrigImages(:,:,:,1);
    plot_Image(handles,currentimage)
    if strcmp(get(handles.tools, 'visible'), 'off');
        set(handles.tools,'visible','on')
    end
elseif status==3 %Convert2binary
    set(handles.fileselector, 'value',1)
    currentimage=handles.bw2(:,:,:,1);
    plot_Image(handles,currentimage)
    if strcmp(get(handles.tools, 'visible'), 'off');
        set(handles.tools,'visible','on')
    end
elseif status==4%Binary_all_Frames
    set(handles.fileselector, 'value',1)
    currentimage=handles.bw2(:,:,:,1);
    plot_Image(handles,currentimage)
    
elseif status==5%Detect_Larvae_in_all_frames
    set(handles.fileselector, 'value',1)
    currentimage=handles.bw2(:,:,:,1);
    plot_Image(handles,currentimage)
    if strcmp(get(handles.tools, 'visible'), 'off');
        set(handles.tools,'visible','on')
    end
    hold on;plot(handles.centroid(1, 1).x,handles.centroid(1, 1).y,'ro');hold off
    if strcmp(get(handles.tools, 'visible'), 'off');
        set(handles.tools,'visible','on')
    end
    hold off
elseif status==6 %Track_Larvae_Paths
    set(handles.tools,'visible','off')
    n_tracks = numel(handles.tracks);
    colors =lines(n_tracks);
    for i_track = 1 : n_tracks
        
        % We use the adjacency tracks to retrieve the points coordinates. It
        % saves us a loop.
        all_points=vertcat(handles.points{:});
        track = handles.adjacency_tracks{i_track};
        track_points = all_points(track, :);
        larvae_tracks{i_track}=[track_points(:,1)  track_points(:, 2)];
        plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :),'linewidth',2)
        hold all
    end
    time_consecutive_larvae=handles.time_consecutive_larvae;
    
    CM = hsv(handles.nFrames);
    for j=1:handles.nFrames
        plot(handles.points{j, 1}(:,1),handles.points{j, 1}(:,end),'MarkerFaceColor',CM(j,:),'MarkerEdgeColor','k','Marker','o','LineStyle','none');
        hold on
        set(gca,'YDir','reverse')
    end
    
    xlabel('X, in pixels')
    ylabel('Y, in pixels')
    hold off
elseif status==7%Calibration
    set(handles.fileselector, 'value',1)
    currentimage=handles.OrigImages(:,:,:,1);
    plot_Image(handles,currentimage)
    set(handles.tools,'visible','on')
    set(handles.Pixel2cm_ratio_text,'string',['The pixel by cm ratio is equal to ' num2str(round((handles.Pix_Cmx)*10)/10) ' and '  num2str(round((handles.Pix_Cmy)*10)/10) ' in x and y, respectively.' ] )
elseif status==8
    set(handles.tools,'visible','off')
end
cla(handles.picture,'reset')
set(handles.picture,'Visible','off');
switchui(['multip0',num2str(status)])

%     sliderdisp(hObject,eventdata,handles)



end

% --- Executes on button press in ExportResults_button.
function ExportResults_button_Callback(hObject, eventdata, handles)
mean_vel_x_and_y=handles.mean_vel_x_and_y;
%mean_Vel_mag=sqrt((mean_vel_x_and_y(:,1).^2)+(mean_vel_x_and_y(:,2).^2));
mean_Vel_mag=handles.mean_Vel_mag;
mean_Distance_x_y=handles.mean_Distance_x_y;
mean_time_span=handles.mean_time_span;
%%
Larvae_velocity_x_and_y=handles.Larvae_velocity_x_and_y;
Distance_x_y=handles.Distance_x_y;
time_span=handles.time_span;
%% Mean results
Mean_col_header={'Larvae_ID','Mean Vx (cm/s)','Mean Vy(cm/s)','Mean Velocity magnitude (cm/s)','Mean distance in x (cm)','Mean distance in y (cm)','Mean time span (seconds)'};%Row cell array (for column labels)
larvaeId=(1:1:length(mean_Vel_mag))';
Data=num2cell([larvaeId mean_vel_x_and_y mean_Vel_mag mean_Distance_x_y mean_time_span]);
Mean_output_matrix=[Mean_col_header;Data];

Larvae_velocity_x_and_y=handles.Larvae_velocity_x_and_y;
Distance_x_y=handles.Distance_x_y;
time_span=handles.time_span;
%col_header={'Larvae_ID','Vx (cm/s)','Vy+(cm/s)','Vy+(cm/s)','Vel Mag(cm/s)'};
col_header={'Larvae_ID','Vx (cm/s)','Vy(cm/s)','Distance in x (cm)','Distance in y (cm)','Average time span (seconds)'};%Row cell array (for column labels)
larvaeId=[];
for i=1:length(handles.Larvae_velocity_x_and_y)
    larvaeId=[larvaeId;repmat(i,size(Larvae_velocity_x_and_y{i,1},1),1)];
end
Data=num2cell([larvaeId cell2mat(Larvae_velocity_x_and_y) cell2mat(Distance_x_y) cell2mat(time_span)]);
output_matrix=[col_header; Data];     %Join cell arrays
xlswrite([handles.filepath,'Results_Larvae_Swim.xls'],output_matrix,'Step_by_Step');     %Write data and both headers
xlswrite([handles.filepath,'Results_Larvae_Swim.xls'],Mean_output_matrix,'Mean_movement_by_Larvae');     %Write data and both headers
Message= msgbox('Done','Done','help');
pause(0.8)
if exist('Message', 'var')
    delete(Message);
    clear('Message');
end
end


% --------------------------------------------------------------------
function Save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,filetype] = uiputfile({'*.tif';'*.jpg';'*.pdf';'*.eps';'*.fig'},'Save plot or image as');
saveDataName = fullfile(PathName,FileName);
fh = figure;
copyobj(handles.picture, fh);
%axes(handles.picture);
saveas(fh, saveDataName); close(fh);
end

% --- Executes on button press in Delete_larvae_in_area.
function Delete_larvae_in_area_Callback(hObject, eventdata, handles)
%%==========
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First the user selects a region, larvae in this region will be deleted.
done=false;
while ~done
    h = msgbox('Please select area containing larvae to be deleted');
    [subIm,roi]=imcrop(handles.picture);
    %wait for user to select region of interest (click-drag rubberband box)
    %determine corners of selected region
    x1=roi(1);  x2=roi(1)+roi(3);  xx=[x1 x1 x2 x2 x1];
    y1=roi(2);  y2=roi(2)+roi(4);  yy=[y1 y2 y2 y1 y1];
    %draw rectangle to indicate selected region
    hLine=line(xx,yy);
    title('Click inside to keep, outside to select again')
    pt=ginput(1);
    %wait for user to keep or discard selected region
    % check: is point inside region?
    if pt(1)>=roi(1) && pt(1)<=roi(1)+roi(3) && pt(2)>=roi(2) && pt(2)<=roi(2)+roi(4)
        done=true;  %exit loop
    else
        delete(hLine)   %remove previous selection
        title('Click-drag-release to select region of interest')
        %loop around to select another region
    end
end
title('The area was detected, larvae located in this area will be deleted');delete(h)
%%
centroid=handles.centroid(selected, 1);
%%delete centroid handles for selected picture frame
handles.centroid(selected, 1).x=[];
handles.centroid(selected, 1).y=[];
Particles_to_Remove=find(centroid.x>x1&centroid.x<x2&centroid.y>y1&centroid.y<y2);
if ~isempty(Particles_to_Remove)
    %Remove selected larvae
    centroid.x(Particles_to_Remove)=[];
    centroid.y(Particles_to_Remove)=[];
    %Re-write handdles
    handles.centroid(selected, 1).x=centroid.x;
    handles.centroid(selected, 1).y=centroid.y;
    guidata(hObject, handles);% Update handles structure
    popupmenu_BackgroundImage_Callback(hObject, eventdata, handles) %Re-plot
end

end

% --------------------------------------------------------------------
function ValidatePaths_menu_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function popupmenu_BackgroundImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_BackgroundImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function dt_multip04_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt_multip04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function dt_multip04_Callback(hObject, eventdata, handles)
% hObject    handle to dt_multip04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_multip04 as text
%        str2double(get(hObject,'String')) returns contents of dt_multip04 as a double
end



function max_linking_distance_Callback(hObject, eventdata, handles)
% hObject    handle to max_linking_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_linking_distance as text
%        str2double(get(hObject,'String')) returns contents of max_linking_distance as a double
end

% --- Executes during object creation, after setting all properties.
function max_linking_distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_linking_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function max_gap_closing_Callback(hObject, eventdata, handles)
% hObject    handle to max_gap_closing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_gap_closing as text
%        str2double(get(hObject,'String')) returns contents of max_gap_closing as a double
end

% --- Executes during object creation, after setting all properties.
function max_gap_closing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_gap_closing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Num_larv_Callback(hObject, eventdata, handles)
% hObject    handle to Num_larv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Num_larv as text
%        str2double(get(hObject,'String')) returns contents of Num_larv as a double
end

% --- Executes during object creation, after setting all properties.
function Num_larv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_larv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in Detect_manually_button.
function Detect_manually_button_Callback(hObject, eventdata, handles)
hold on
%counter=0;
Num_larv=str2double(get(handles.Num_larv,'string')) ;
if Num_larv>0
    for j=1:Num_larv
        [x(j,1),y(j,1)] = ginput(1);
        plot(x(j,1),y(j,1),'bo','MarkerSize',5)
        %pause
    end
else
    return
end
%[x y]=ginput(1);
plot(x,y,'bo')
hold off

%%==========
toggler=retr('toggler');
selected=floor(get(handles.fileselector, 'value'))-toggler*(1);
if selected==0 % If return to initial frame while toggle is activated, setup selected to first frame
    selected=1;
end

handles.centroid(selected, 1).x=vertcat(handles.centroid(selected, 1).x,x);
handles.centroid(selected, 1).y=vertcat(handles.centroid(selected, 1).y,y);
guidata(hObject, handles);% Update handles structure
hold off
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double
end

% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double
end

% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double
end

% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
cla(handles.picture,'reset');
mean_Vel_mag=handles.mean_Vel_mag;
%Larvae_velocity_x_and_y=cell2mat(handles.Larvae_velocity_x_and_y);
mean_vel_x_and_y=handles.mean_vel_x_and_y;
guidata(hObject, handles);% Update handles structure

if (get(handles.popupmenu_Display,'value'))==2
    cla(handles.picture,'reset');
    figure('color','w')
end
set(handles.picture,'Visible','on');
if (get(handles.popupmenu_plot_type,'value'))==1
    if (get(handles.Popup_variable,'value'))==1
        bar(mean_Vel_mag,'b');
        xlabel('Larvae')
        ylabel('Velocity magnitude [cm/s]')
    elseif (get(handles.Popup_variable,'value'))==2
        barh(mean_vel_x_and_y(:,1),'g');
        ylabel('Larvae')
        xlabel('Horizontal component of velocity [cm/s]')
    elseif (get(handles.Popup_variable,'value'))==3
        bar(mean_vel_x_and_y(:,2),'c');
        xlabel('Larvae')
        ylabel('Vertical component of velocity [cm/s]')
    end
end
%
if (get(handles.popupmenu_plot_type,'value'))==2
    if (get(handles.Popup_variable,'value'))==1
        [counts, binValues] = hist(mean_Vel_mag);
        bar(binValues, counts/sum(counts),'b', 'barwidth', 1);
        xlabel('Velocity magnitude [cm/s]')
        ylabel('Percentage of larvae')
    elseif (get(handles.Popup_variable,'value'))==2
        [counts, binValues] = hist(mean_vel_x_and_y(:,1));
        bar(binValues, counts/sum(counts),'g', 'barwidth', 1);
        xlabel('Horizontal component of velocity [cm/s]')
        ylabel('Percentage of larvae')
    elseif (get(handles.Popup_variable,'value'))==3
        [counts, binValues] = hist(mean_vel_x_and_y(:,2));
        bar(binValues, counts/sum(counts),'c', 'barwidth', 1);
        xlabel('Vertical component of velocity [cm/s]')
        ylabel('Percentage of larvae')
    end
end %if hystogram
%
if (get(handles.popupmenu_plot_type,'value'))==3
    if (get(handles.Popup_variable,'value'))==1
        scatter([1:1:length(mean_Vel_mag(:,1))],mean_Vel_mag(:,1),150,'b','fill');
        xlabel('Larvae')
        ylabel('Velocity magnitude [cm/s]')
    elseif (get(handles.Popup_variable,'value'))==2
        scatter([1:1:length(mean_vel_x_and_y(:,1))],mean_vel_x_and_y(:,1),150,'g','fill');
        xlabel('Larvae')
        ylabel('Horizontal component of velocity [cm/s]')
    elseif (get(handles.Popup_variable,'value'))==3
        scatter([1:1:length(mean_vel_x_and_y(:,2))],mean_vel_x_and_y(:,2),150,'c','fill');
        xlabel('Larvae')
        ylabel('Vertical component of velocity [cm/s]')
    end
end



set(handles.text_Stats,'String',['Summary stats: ','Mean velocity magnitude=',num2str(round(mean(mean_Vel_mag)*10)/10),' cm/s'])
end


function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double
end

% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
end

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit35_Callback(hObject, ~, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double
end

% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu_plot_type.
function popupmenu_plot_type_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function popupmenu_plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu_Display.
function popupmenu_Display_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Display
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
end


% --------------------------------------------------------------------
function Save_menu_Callback(hObject, eventdata, handles)
if isfield(handles,'filepath')
    filepath=handles.filepath;
else
    errordlg('Please start a new sesion','Error');
    return
end
% if isempty(filepath)
%     sessionpath=retr('pathname');
% end
[Session_FileName,Session_PathName] = uiputfile('*.mat','Save current session as...',fullfile(filepath,'LarvaeSwim.mat'));
% if isequal(FileName,0) | isequal(PathName,0)
% else
%     put('sessionpath',PathName );
handles.Session_FileName=Session_FileName;
handles.Session_PathName=Session_PathName;
savesessionfuntion(handles)
% end
guidata(hObject, handles);% Update handles structure
end


function savesessionfuntion (handles)
ed=msgbox('Please wait, saving session...','Saving','help');
Session_FileName=handles.Session_FileName;
Session_PathName=handles.Session_PathName;

%Newer versions of Matlab do really funny things when the following vars are not empty...:
app.GUIDEOptions =[];
app.GUIOnScreen  =[];
app.Listeners  =[];
app.SavedVisible  =[];
app.ScribePloteditEnable  =[];
app.UsedByGUIData_m  =[];
app.lastValidTag =[];
iptPointerManager=[];
try
    iptPointerManager(gcf, 'disable')
end

save('-v7.3', fullfile(Session_PathName,Session_FileName), '-struct', 'app')
clear app hgui  iptPointerManager

%% From GUI
try
    radius=str2double(get(handles.radius,'String'));
    threshold=str2double(get(handles.threshold,'String'));
    P=str2double(get(handles.P,'String'));
    gaps=str2double(get(handles.gaps,'String'));
    
    %% From Handles
    filepath=handles.filepath;
    filename=handles.filename;
    nFrames=handles.nFrames;
    status=handles.status;
    OrigImages=handles.OrigImages;
    Crop_OrigImages=handles.Crop_OrigImages;
    Images=handles.Images;
    bw2=handles.bw2;
    subIm=handles.subIm;
    centroid=handles.centroid;
    roi=handles.roi;
    a=handles.a;
    larvae_tracks=handles.larvae_tracks;
    time_consecutive_larvae=handles.time_consecutive_larvae;
    points=handles.points;
    track_points=handles.track_points;
    tracks=handles.tracks;
    adjacency_tracks=handles.adjacency_tracks;
    Pix_Cmx=handles.Pix_Cmx;
    Pix_Cmy=handles.Pix_Cmy;
    mean_Vel_mag=handles.mean_Vel_mag;
    Larvae_velocity_x_and_y=handles.Larvae_velocity_x_and_y;
    mean_vel_x_and_y=handles.mean_vel_x_and_y;
    
    
end
clear handles
if exist('ed', 'var')
    delete(ed);
    clear('ed');
end
save('-v7.3', fullfile(Session_PathName,Session_FileName),'-append');
ed=msgbox('Done','Saving','help');
if exist('ed', 'var')
    delete(ed);
    clear('ed');
end
end
function dt_multip05_Callback(hObject, eventdata, handles)
end
function dt_multip05_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function Settings_Callback(hObject, eventdata, handles)

end
function Larvae_Tracking_Menu_Callback(hObject, eventdata, handles)

end
function Larvae_Identificat_Callback(hObject, eventdata, handles)
% hObject    handle to Larvae_Identificat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
function filenamebox_Callback(hObject, eventdata, handles)
% hObject    handle to filenamebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filenamebox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filenamebox
end
function filenamebox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function loadimgbutton_Callback(hObject, eventdata, handles)
%
end
function pushbutton12_Callback(hObject, eventdata, handles)

end
function edit9_Callback(hObject, eventdata, handles)
end
function edit9_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit12_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit13_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit14_Callback(hObject, eventdata, handles)
end

function edit14_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%useful code for future

% % maximum distance from line to the mouse pointer (arbitrary)
%     DISTANCE_THRESHOLD = 2;
%     % get the handles structure
%     %handles = guidata(figHandle);
% %     % get the position where the mouse button was pressed (not released)
%     % within the GUI
%     [x,y] = ginput(1)
%     currentPoint = get(handles.picture, 'CurrentPoint');
%     x            = currentPoint(1,1);
%     y            = currentPoint(1,2);
%     % get the position of the axes within the GUI
%     axesPos = get(handles.axes1,'Position');
%     minx    = axesPos(1);
%     miny    = axesPos(2);
%     maxx    = minx + axesPos(3);
%     maxy    = miny + axesPos(4);
%     % is the mouse down event within the axes?
%     if x>=minx && x<=maxx && y>=miny && y<=maxy
%         % do we have graphics objects?
%         if isfield(handles,'plotHandles')
%             % get the position of the mouse down event within the axes
%             currentPoint = get(handles.axes1, 'CurrentPoint');
%             x            = currentPoint(2,1);
%             y            = currentPoint(2,2);
%             % we are going to use the x and y data for each graphic object
%             % and determine which one is closest to the mouse down event
%             minDist      = Inf;
%             minHndlIdx   = 0;
%            for k=1:length(handles.plotHandles)
%                xData = get(handles.plotHandles(k),'XData');
%                yData = get(handles.plotHandles(k),'YData');
%                dist  = min((xData-x).^2+(yData-y).^2);
%                if dist<minDist && dist<DISTANCE_THRESHOLD
%                    minHndlIdx = k;
%                    minDist = dist;
%                end
%            end
%            % if we have a graphics handle that is close to the mouse down
%            % event/position, then save its index into the plotHandles array
%            if minHndlIdx~=0
%                handles.graphicToDeleteHandleIdx = minHndlIdx;
%            else
%                handles.graphicToDeleteHandleIdx = [];
%            end
%            % change the line style of the selected object
%            for k=1:length(handles.plotHandles)
%                if k==minHndlIdx
%                    set(handles.plotHandles(k),'LineStyle',':');
%                else
%                    set(handles.plotHandles(k),'LineStyle','-');
%                end
%            end
%            guidata(figHandle,handles);
%         end
%     end
% end

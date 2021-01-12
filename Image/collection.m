function varargout = collection(varargin)
% COLLECTION MATLAB code for collection.fig
%      COLLECTION, by itself, creates a new COLLECTION or raises the existing
%      singleton*.
%
%      H = COLLECTION returns the handle to a new COLLECTION or the handle to
%      the existing singleton*.
%
%      COLLECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLLECTION.M with the given input arguments.
%
%      COLLECTION('Property','Value',...) creates a new COLLECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before collection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to collection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help collection

% Last Modified by GUIDE v2.5 08-May-2016 13:29:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @collection_OpeningFcn, ...
                   'gui_OutputFcn',  @collection_OutputFcn, ...
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


% --- Executes just before collection is made visible.
function collection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to collection (see VARARGIN)

% Choose default command line output for collection
handles.output = hObject;

% add path
addpath('utilities', 'group_transfer');
warning('off','all');
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes collection wait for user response (see UIRESUME)
% uiwait(handles.figure1);
        
% --- Outputs from this function are returned to the command line.
function varargout = collection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'units','normalized');


% Initialize parameters
handles.n = 0;
handles.CollectionGridSize = [2 6];
handles.Margin = [20 20];
handles.HighlightBrochure = 0;
handles.HighlightRef = 0;
handles.HighlightOutput = 0;
handles.mode = 4;
handles.ExtractPalette = 0;
handles.NotConnectedIcon = im2double(imread('icon\not_connect.png'));
handles.ConnectedIcon = im2double(imread('icon\connect.png'));
handles.LockIcon = im2double(imread('icon\lock.png'));
handles.BackgroundColor = 0.5; 
handles.K = 5; % 5 is standard using for drawing palette

axesPos = getpixelposition(handles.axesImageCollection);
W = uint16(axesPos(3)); H = uint16(axesPos(4));
handles.CollectionImageSize = [H, W];
subH = uint16((H - (handles.CollectionGridSize(1)-1)*handles.Margin(1))/handles.CollectionGridSize(1));
subW = uint16((W - (handles.CollectionGridSize(2)-1)*handles.Margin(2))/handles.CollectionGridSize(2));
handles.subImageSize = [subH, subW];

handles.BrochurePalette = [];
handles.ConstraintFlag = 0;
handles.VotingMode = 0;
handles.RecolorBrochure = 0;
handles.sigma = 80;
handles.Pref = [];
handles.isNewPalette = 1;
load('data.mat');
handles.LUT1 = LUT1;
handles.LUT2 = LUT2;

% parameters for optimization
handles.AutomaticTransfer = 0;
handles.solver_no_group = 3;
handles.solver_gamma = 0;
handles.solver_eta = 10^10;

% Display Brochure
% displayBrochureImage(hObject, handles);

% Output image
axesPos = getpixelposition(handles.axesOutputImage);
W = floor(axesPos(3)); H = floor(axesPos(4));
handles.Iout = handles.BackgroundColor *ones(H, W, 3);
axes(handles.axesOutputImage);
h = imshow(handles.Iout);

% Reference palette
axesPos = getpixelposition(handles.axesPalette);
W = floor(axesPos(3)); H = floor(axesPos(4));
handles.PaletteSize = [H, W];
PIref = handles.BackgroundColor *ones(H, W, 3);
axes(handles.axesPalette);
imshow(PIref);

% Image Collection
displayImageCollection(hObject, handles);

% Update handles structure
guidata(hObject, handles);

statusbar(gcf, 'Done'); 
% set(get(gcf,'JavaFrame'),'Maximized',1);


% --------------------------------------------------------------------
function displayBrochureImage(hObject, handles)
axesPos = getpixelposition(handles.axesBrochureImage);
W = floor(axesPos(3)); H = floor(axesPos(4));
rs = W/H;
Ishow = handles.BackgroundColor *ones(H, W, 3);

I = handles.BrochureImageDisplay;
ri = size(I,2)/size(I,1);
if (ri > rs)
    Wi = W; Hi = floor(Wi/ri);
    iw = 0; ih = floor((H-Hi)/2);
else
    Hi = H;
    Wi = floor(Hi*ri);
    iw = floor((W-Wi)/2); ih = 0;
end
I = imresize(I, [Hi, Wi]);
Ishow(ih+1:ih+Hi, iw+1:iw+Wi, :) = I;
axes(handles.axesBrochureImage);
imshow(Ishow);
    
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function handles = displayImageCollection(hObject, handles)
% Create image for the collection
m = handles.CollectionGridSize(1); n = handles.CollectionGridSize(2);
axesPos = getpixelposition(handles.axesImageCollection);
W = uint16(axesPos(3)); H = uint16(axesPos(4));
outSize = 2;

handles.CollectionImageSize = [H, W];
Hsub = uint16((H - 4*outSize - (handles.CollectionGridSize(1)-1)*handles.Margin(1))/m);
Wsub = uint16((W - 4*outSize - (handles.CollectionGridSize(2)-1)*handles.Margin(2))/n);


Hmar = handles.Margin(1); Wmar = handles.Margin(2);
H = handles.CollectionImageSize(1); W = handles.CollectionImageSize(2);

rs = double(Wsub)/double(Hsub);

Icol = zeros(H, W, 3);
c = get(handles.figure1, 'Color');
for i = 1:3
    Icol(:,:,i) = c(i);
end

if handles.AutomaticTransfer == 1
    unmatch = checkUnMatchBrochurePalette(handles.PrefMatch, handles.Match);
end
% Get color of main figure
% c = get(handles.figure1, 'Color');
c = handles.BackgroundColor *ones(1,3);

outColor = [1, 1, 1];
for i = 1:m
    for j = 1:n
        id = (i-1)*n + j;
        hs = outSize + Hmar*(i-1) + Hsub*(i-1);
        ws = outSize + Wmar*(j-1) + Wsub*(j-1);
        
        % insert highlight for output image
        if(id <= handles.n && id == handles.Idout)
            for b = 1:3
                Icol(hs+1-outSize:hs+Hsub+outSize, ws+1-outSize:ws+Wsub+outSize, b) = outColor(b);
            end
        end

        % insert background 
        for b = 1:3
            Icol(hs+1:hs+Hsub, ws+1:ws+Wsub, b) = c(b);
        end
        
        % insert image       
        if (id > handles.n)
            continue;
        end

        if get(handles.ckbShowOriginal, 'Value') == 0
            I = handles.Imod{id};
        else
            I = handles.Iorg{id};
        end
        
        ri = size(I,2)/size(I,1);
        if (ri > rs)
            Wi = Wsub; Hi = floor(Wi/ri);
            iw = 0; ih = floor((Hsub-Hi)/2);
        else
            Hi = Hsub;
            Wi = floor(Hi*ri);
            iw = floor((Wsub-Wi)/2); ih = 0;
        end
        I = imresize(I, [Hi, Wi]);
        Icol(hs+ih+1:hs+ih+Hi, ws+iw+1:ws+iw+Wi, :) = I;
        
        if(handles.Reference(id))
            [iconh, iconw, ~] = size(handles.LockIcon);
            Icol(hs+Hsub-iconh+1:hs+Hsub, ws+Wsub-iconw+1:ws+Wsub, :) = handles.LockIcon;
        end
    end
end
axes(handles.axesImageCollection);
h = imshow(Icol);
set(h, 'ButtonDownFcn',@subImageCollection_Callback)




function subImageCollection_Callback(o,e)
%# show selected image in a new figure
handles = guidata(o);

curPoint = get(handles.axesImageCollection, 'CurrentPoint');
step = double(handles.CollectionImageSize)./double(handles.CollectionGridSize);
m = ceil(curPoint(1,2)/step(1));
n = ceil(curPoint(1,1)/step(2));
idx = (m-1)*handles.CollectionGridSize(2)+n;
if idx > handles.n
    return;
end

mode = get(gcf,'Selectiontype');
if strcmp(mode, 'alt')
    handles.Reference(idx) = 1 - handles.Reference(idx);
    
    n = 0;
    handles.W3all = zeros(256,1);
    handles.H3all = zeros(256,2);
    for i = 1:handles.n
        if(handles.Reference(i))
            handles.W3all = handles.W3all + handles.W3{i};
            handles.H3all = handles.H3all + repmat(handles.W3{i}, 1, 2).*handles.H3{i};
            n = n + 1;
        end
    end
    
    if (n == 0)  
        for i = 1:handles.n
            handles.W3all = handles.W3all + handles.W3{i};
            handles.H3all = handles.H3all + repmat(handles.W3{i}, 1, 2).*handles.H3{i};
        end               
    end
    % Normalize total histogram
    map = handles.W3all~=0;
    handles.W3all = handles.W3all(map);
    handles.H3all = handles.H3all(map,:)./repmat(handles.W3all,1,2);
    [handles.Pmod, handles.Match, handles.Pref] = solve_optimal_all_palette(handles.Porg, handles.Reference, handles.Horg, handles.H3all, handles.W3all, handles.solver_no_group, handles.solver_gamma, handles.solver_eta);
    handles.PrefMatch = zeros(length(handles.Pref), 1);
    
    % Transfer color
    if handles.AutomaticTransfer
        statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
        for i = 1:handles.n
            [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
            statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
        end
        displayOutputImage(handles);
		statusbar('Done');
    end
end
    
% display palette
handles.Idout = idx;
displayOutputImage(handles);    

handles.isNewPalette = 1;

handles = displayPalette(handles.figure1, handles);
% Update image collection
handles = displayImageCollection(handles.figure1, handles);

guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function handles = displayPalette(hObject, handles)

axesCur = handles.axesPalette;
Pbro = handles.BrochurePalette;


if handles.ExtractPalette
    Pref = handles.Pref;
    Pout = handles.Porg{handles.Idout};
else 
    Pref = [];
    Pout = [];
end

if handles.isNewPalette == 1
    handles.HighlightBrochure = 0;
    handles.HighlightRef = 0;
    handles.HighlightOutput = 0;
end


% Create palette image 
margin = 10;
Hsub = 30;
Wsub = 45;
mref = size(Pref, 1);
mout = size(Pout, 1);
mbro = size(Pbro, 1);
mcell = max([mref, mout, mbro,5]);
H = mcell*Hsub + (mcell+1)*margin;

axesPos = getpixelposition(axesCur);
rs = axesPos(3)/axesPos(4);
W = round(H * rs);
Pimage = handles.BackgroundColor*ones(H, W, 3);

% Insert brochure palette
highlightColor = [0 0 1];
norColor = [1 1 1]; 
for i = 1:mbro
    is = (i-1)*(Hsub+margin) + margin;  
    for j = 1:3      
        if handles.HighlightBrochure == i
            Pimage(is+1:is+Hsub, margin+1:Wsub+margin, j) = highlightColor(j);
            Pimage(is+3:is+Hsub-2, margin+3:Wsub-2+margin, j) = Pbro(i,j);  
        else
            Pimage(is+1:is+Hsub, margin+1:Wsub+margin, j) = norColor(j);
            Pimage(is+3:is+Hsub-2, margin+3:Wsub-2+margin, j) = Pbro(i,j);
        end
    end
end

% Insert group pallete
for i = 1:mref
    is = (i-1)*(Hsub+margin) + margin;  
    for j = 1:3        
        js = (W-Wsub)/2;
        if handles.HighlightRef == i
            Pimage(is+1:is+Hsub, js+1:js+Wsub, j) = highlightColor(j);
            Pimage(is+3:is+Hsub-2, js+3:js+Wsub-2, j) = Pref(i,j);
        else
            Pimage(is+1:is+Hsub, js+1:js+Wsub, j) = norColor(j);
            Pimage(is+3:is+Hsub-2, js+3:js+Wsub-2, j) = Pref(i,j);
        end       
    end
end

% Insert current viewed image palette
for i = 1:mout
    for j = 1:3
        is = (i-1)*(Hsub+margin) + margin;         
        
        if handles.HighlightOutput == i
            Pimage(is+1:is+Hsub, W-Wsub+1-margin:W-margin, j) = highlightColor(j);
            Pimage(is+3:is+Hsub-2, W-Wsub+3-margin:W-2-margin, j) = Pout(i,j);
        else
            Pimage(is+1:is+Hsub, W-Wsub+1-margin:W-margin, j) = norColor(j);
            Pimage(is+3:is+Hsub-2, W-Wsub+3-margin:W-2-margin, j) = Pout(i,j);
        end
        
    end
end
    
% show palette
handles.PaletteSize = size(Pimage);
axes(axesCur);
h = imshow(Pimage);
set(h, 'ButtonDownFcn',@subPalette_Callback)
hold on;
lineColor = [1, 1, 1];
% draw all line from Brochure to Reference
for i = 1:mref
    if handles.PrefMatch(i) > 0
        ie = (i-1)*(Hsub+margin)+margin;
        is = (handles.PrefMatch(i)-1)*(Hsub+margin)+margin;
        je = (W - Wsub)/2;
        line([Wsub+6+margin, je-5], [is + Hsub/2, ie+Hsub/2], 'Color', lineColor, 'LineWidth', 4);
    end
end

% draw all line from Reference to individual
for i = 1:mout
    if handles.Match{handles.Idout}(i) > 0
        ie = (i-1)*(Hsub+margin)+margin;
        is = (handles.Match{handles.Idout}(i)-1)*(Hsub+margin)+margin;
        js = (W - Wsub)/2;
        line([js+Wsub+6, W-Wsub-5-margin], [is + Hsub/2, ie+Hsub/2], 'Color', lineColor, 'LineWidth', 4);
    end
end
hold off;


% --- Executes on button press in subPalettes.
function subPalette_Callback(o, e)
handles = guidata(o);

if handles.ExtractPalette == 0
    return;
end

curPoint = get(handles.axesPalette, 'CurrentPoint');
step = handles.PaletteSize(1)/max([size(handles.BrochurePalette,1), size(handles.Pref,1), size(handles.Porg{handles.Idout},1), handles.K]); % 5 is standard
m = ceil(curPoint(1,2)/step);

n = curPoint(1,1);
if (n <= handles.PaletteSize(2)/3)
    m = max(min(m, size(handles.BrochurePalette, 1)), 1);
    if handles.HighlightBrochure == m
        handles.HighlightBrochure = 0;
    else
        handles.HighlightBrochure = m;
    end
elseif (n>= handles.PaletteSize(2)/3 && n <= 2*handles.PaletteSize(2)/3)
    m = max(min(m, size(handles.Pref, 1)), 1);
    mode = get(gcf,'Selectiontype');
    if strcmp(mode, 'alt') 
        c = uisetcolor;
        if(~isequal(c, 0))
            handles.Pref(m,:) = c;            
            handles.Pmod = assignColorReference2InputManual(handles.Porg, handles.Match, handles.Pref);
            % Transfer color
            if handles.AutomaticTransfer
                statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
                for i = 1:handles.n
                    [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
                    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
                end
                displayOutputImage(handles);
                statusbar('Done');
            end
        end
    else
        if handles.HighlightRef == m
            handles.HighlightRef = 0;
        else
            handles.HighlightRef = m;
        end
    end
elseif (n >= 2*handles.PaletteSize(2)/3)
    m = max(min(m, handles.K), 1);
    if handles.HighlightOutput == m
        handles.HighlightOutput = 0;
    else
        handles.HighlightOutput = m;
    end
else
    return;
end

if handles.HighlightBrochure > 0 && handles.HighlightRef > 0
    if handles.PrefMatch(handles.HighlightRef) == handles.HighlightBrochure
        handles.PrefMatch(handles.HighlightRef) = 0;
        handles.Prefbro(handles.HighlightRef,:) = handles.Pref(handles.HighlightRef,:);
    else
        handles.PrefMatch(handles.HighlightRef) = handles.HighlightBrochure;
        handles.Prefbro(handles.HighlightRef,:) = handles.BrochurePalette(handles.HighlightBrochure,:);
    end
    handles.isNewPalette = 1;   
    handles.Pmod = assignColorReference2InputManual(handles.Porg, handles.Match, handles.Prefbro);
    
    % Transfer color
    if handles.AutomaticTransfer
        statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
        for i = 1:handles.n
            [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
            statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
        end
        displayOutputImage(handles);
		statusbar('Done');
    end
    
    
elseif handles.HighlightRef > 0 && handles.HighlightOutput > 0
    % Toggle the match
    if handles.Match{handles.Idout}(handles.HighlightOutput) == handles.HighlightRef
        handles.Match{handles.Idout}(handles.HighlightOutput) = 0;        
        handles.Pmod{handles.Idout}(handles.HighlightOutput, :) = handles.Porg{handles.Idout}(handles.HighlightOutput, :);
    else
        handles.Match{handles.Idout}(handles.HighlightOutput) = handles.HighlightRef;
        handles.Pmod{handles.Idout}(handles.HighlightOutput, :) = handles.Pref(handles.HighlightRef, :);
    end
    handles.isNewPalette = 1;
    idx = handles.Idout;
    if handles.AutomaticTransfer
        [handles.Imod{idx}, ~] = lab_transfer(handles.Iorg{idx}, handles.Porg{idx}, handles.Pmod{idx}, handles.LUT1, handles.LUT2);
    end
    displayOutputImage(handles);    
else 
    if handles.HighlightBrochure == 0 && handles.HighlightRef == 0
        displayOutputImage(handles);
    end
    handles.isNewPalette = 0;
end

handles = displayPalette(handles.figure1, handles);
handles = displayImageCollection(handles.figure1, handles);
guidata(handles.figure1, handles);


function handles = ExtractPalette(hObject, handles)
if handles.n == 0
    return;
end

bin = 15;
statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
handles.W3all = zeros((bin+1)^2, 1);
handles.H3all = zeros((bin+1)^2, 2);
handles.Porg = [];
handles.Horg = [];
for i = 1:handles.n
    disp(['Extracting palette image ' num2str(i)]);
    [handles.Porg{i}, handles.Horg{i}, handles.Labels{i}, handles.W3{i}, handles.H3{i}] = extract_theme_elbow(handles.Iorg{i}, 2, handles.LUT1, bin);
    handles.Match{i} = zeros(length(handles.Porg{i}),1);
    handles.W3all = handles.W3all + handles.W3{i};
    handles.H3all = handles.H3all + repmat(handles.W3{i}, 1, 2).*handles.H3{i};
    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
end
% Normalize total histogram
map = handles.W3all~=0;
handles.W3all = handles.W3all(map);
handles.H3all = handles.H3all(map,:)./repmat(handles.W3all,1,2);
% Copy the reference palette
handles.isNewPalette = 1;
handles.ExtractPalette = 1;
% Update handles structure
handles = displayPalette(hObject, handles);
statusbar(gcf, 'Done');




function displayOutputImage(handles)
if get(handles.ckbShowOriginal, 'Value') == 1
    I = handles.Iorg{handles.Idout};
else
    I = handles.Imod{handles.Idout};
end
if handles.HighlightRef == 0 && handles.HighlightOutput > 0
    labels = handles.Labels{handles.Idout};
    mask = labels == handles.HighlightOutput;
    I = I.*repmat(mask, [1 1 3]);
end



axesPos = getpixelposition(handles.axesOutputImage);
W = floor(axesPos(3)); H = floor(axesPos(4));
rs = W/H;
Ishow = handles.BackgroundColor *ones(H, W, 3);
% c = get(handles.figure1, 'Color');
% for i = 1:3
%     Ishow(:,:,i) = c(i);
% end

ri = size(I,2)/size(I,1);
if (ri > rs)
    Wi = W; Hi = floor(Wi/ri);
    iw = 0; ih = floor((H-Hi)/2);
else
    Hi = H;
    Wi = floor(Hi*ri);
    iw = floor((W-Wi)/2); ih = 0;
end
I = imresize(I, [Hi, Wi]);
Ishow(ih+1:ih+Hi, iw+1:iw+Wi, :) = I;

axes(handles.axesOutputImage);

imshow(Ishow);






% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:handles.n
    imwrite(handles.Imod{i}, [handles.pathname sprintf('output_%02d.png', i)]);
end
disp('All output files are saved in the same folder with input files.');



% --- Executes on slider movement.
function sliderDistance_Callback(hObject, eventdata, handles)
% hObject    handle to sliderDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if(handles.n <=0)
    return;
end

% Update number of group theme colors
val = round(get(handles.sliderDistance, 'Value'));
set(handles.sliderDistance, 'Value', val);
handles.solver_no_group = 7 - val + 1;
set(handles.txtNumber, 'String', num2str(handles.solver_no_group));

% Compute the new group them colors and apply recoloring
[handles.Pmod, handles.Match, handles.Pref] = solve_optimal_all_palette(handles.Porg, handles.Reference, handles.Horg, handles.H3all, handles.W3all, handles.solver_no_group, handles.solver_gamma, handles.solver_eta);
handles.PrefMatch = zeros(length(handles.Pref), 1);
statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
for i = 1:handles.n
    [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
end

displayOutputImage(handles);
% Display Output palette
handles = displayPalette(hObject, handles);
% Update the image collection
handles = displayImageCollection(hObject, handles);
% Update guidata
guidata(hObject, handles);
statusbar('Done');

% --- Executes during object creation, after setting all properties.
function sliderDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function menuImage_Callback(hObject, eventdata, handles)
% hObject    handle to menuImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.png;*.jpg;*.bmp'},'Select an Image file', '.\images\', 'MultiSelect', 'on');
if(~isequal(filename,0))
    handles.n = length(filename);      
    handles.pathname = pathname;
    for i = 1:handles.n
        handles.Iorg{i} = im2double(imread([pathname filename{i}]));
        handles.Porg{i} = ones(1, 3);
        handles.Imod{i} = handles.Iorg{i};
        handles.Pmod{i} = handles.Porg{i};
        handles.Match{i} = 1;
    end
    
    handles.Pref = [];
    handles.Prefbro = [];
    handles.Reference = zeros(handles.n, 1);
    handles.isNewPalette = 0;
    handles.ExtractPalette = 0;
    
    % Set the number of group theme colors
    handles.solver_no_group = 3;
    set(handles.sliderDistance, 'Value', 7 - handles.solver_no_group + 1);
    set(handles.txtNumber, 'String', num2str(handles.solver_no_group));
    handles.AutomaticTransfer = 0;
    set(handles.chbTransferColor, 'Value', 0);
    
    handles.flag_unmatch = 0;
    handles.flag_collapse = 0;
    set(handles.chb_allow_unmatch, 'Value', 0);
    set(handles.chb_avoid_collapse, 'Value' ,0);
    
    % Default reference
    handles.Idout = 1;

    handles = ExtractPalette(hObject, handles);
    % Default output
    displayOutputImage(handles);   

    % Show image collection
    displayImageCollection(hObject, handles);

    % Update handles structure
    handles.FileType = 2;
    guidata(hObject, handles);

end



% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
% hObject    handle to menuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:handles.n
    imwrite(handles.Imod{i}, [handles.pathname sprintf('output_%02d.png', i)]);
end


% --------------------------------------------------------------------
function menuExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuTools_Callback(hObject, eventdata, handles)
% hObject    handle to menuTools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuImport_Callback(hObject, eventdata, handles)
% hObject    handle to menuImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[P, M, select] = displayBrochureTheme();
% set(handles.sliderRotatePalette, 'Value', 0);
if select > 0
    handles.BrochurePalette = P; 
    handles.BrochurePaletteOrg = P;
    [handles.Prefbro, handles.PrefMatch] = assignColorBrochure2Reference(handles.BrochurePalette, handles.Pref, handles.ConstraintFlag);
    handles.Pmod = assignColorReference2InputManual(handles.Porg, handles.Match, handles.Prefbro);
    handles = displayPalette(hObject, handles);
    handles = displayImageCollection(hObject, handles);
    guidata(hObject, handles); 
end


% --------------------------------------------------------------------
function menuTransfer_Callback(hObject, eventdata, handles)
% hObject    handle to menuTransfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Pmod = assignColorReference2InputManual(handles.Porg, handles.Match, handles.Prefbro);
statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
for i = 1:handles.n
    [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
end

displayOutputImage(handles);
% Display Output palette
handles = displayPalette(hObject, handles);
% Update the image collection
handles = displayImageCollection(hObject, handles);
% Update guidata
guidata(hObject, handles);
statusbar('Done');



% --- Executes on button press in ckbShowOriginal.
function ckbShowOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to ckbShowOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ckbShowOriginal
if handles.n == 0
    return;
end
% Update the image collection
displayImageCollection(hObject, handles);
displayOutputImage(handles);
guidata(hObject, handles);




% --- Executes on button press in chbTransferColor.
function chbTransferColor_Callback(hObject, eventdata, handles)
% hObject    handle to chbTransferColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbTransferColor

if(get(hObject, 'Value') == 1)
    handles.AutomaticTransfer = 1;
    [handles.Pmod, handles.Match, handles.Pref] = solve_optimal_all_palette(handles.Porg, handles.Reference, handles.Horg, handles.H3all, handles.W3all, handles.solver_no_group, handles.solver_gamma, handles.solver_eta);
    handles.PrefMatch = zeros(length(handles.Pref), 1);
    statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
    for i = 1:handles.n
        [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
        statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
    end

    displayOutputImage(handles);
    % Display Output palette
    handles = displayPalette(hObject, handles);
    % Update the image collection
    handles = displayImageCollection(hObject, handles);
    % Update guidata
    guidata(hObject, handles);
    statusbar('Done');

else
    handles.AutomaticTransfer = 0;
end
guidata(hObject, handles);


% --- Executes on button press in chb_avoid_collapse.
function chb_avoid_collapse_Callback(hObject, eventdata, handles)
% hObject    handle to chb_avoid_collapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chb_avoid_collapse
handles.flag_collapse = get(hObject,'Value');
if handles.flag_collapse
    handles.solver_gamma = 100;
    if handles.flag_unmatch
        handles.solver_eta = 25000;
    else
        handles.solver_eta = 10^10;
    end
else
    handles.solver_gamma = 0;
    if handles.flag_unmatch
        handles.solver_eta = 1000;
    else        
        handles.solver_eta = 10^10;
    end
end
    
[handles.Pmod, handles.Match, handles.Pref] = solve_optimal_all_palette(handles.Porg, handles.Reference, handles.Horg, handles.H3all, handles.W3all, handles.solver_no_group, handles.solver_gamma, handles.solver_eta);
handles.PrefMatch = zeros(length(handles.Pref), 1);
statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
for i = 1:handles.n
    [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
end

displayOutputImage(handles);
% Display Output palette
handles = displayPalette(hObject, handles);
% Update the image collection
handles = displayImageCollection(hObject, handles);
% Update guidata
guidata(hObject, handles);
statusbar('Done');

% --- Executes on button press in chb_allow_unmatch.
function chb_allow_unmatch_Callback(hObject, eventdata, handles)
% hObject    handle to chb_allow_unmatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chb_allow_unmatch
handles.flag_unmatch = get(hObject,'Value');
if handles.flag_collapse
    handles.solver_gamma = 100;
    if handles.flag_unmatch
        handles.solver_eta = 25000;
    else
        handles.solver_eta = 10^10;
    end
else
    handles.solver_gamma = 0;
    if handles.flag_unmatch
        handles.solver_eta = 1000;
    else        
        handles.solver_eta = 10^10;
    end
end
    
[handles.Pmod, handles.Match, handles.Pref] = solve_optimal_all_palette(handles.Porg, handles.Reference, handles.Horg, handles.H3all, handles.W3all, handles.solver_no_group, handles.solver_gamma, handles.solver_eta);
handles.PrefMatch = zeros(length(handles.Pref), 1);
statusbar('Processing %d of %d (%.1f%%)...',0,handles.n,0);
for i = 1:handles.n
    [handles.Imod{i}, ~] = lab_transfer(handles.Iorg{i}, handles.Porg{i}, handles.Pmod{i}, handles.LUT1, handles.LUT2);
    statusbar('Processing %d of %d (%.1f%%)...',i,handles.n,100*i/handles.n); 
end

displayOutputImage(handles);
% Display Output palette
handles = displayPalette(hObject, handles);
% Update the image collection
handles = displayImageCollection(hObject, handles);
% Update guidata
guidata(hObject, handles);
statusbar('Done');


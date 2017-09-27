
function varargout = iFextA(varargin)
% IFEXTA MATLAB code for iFextA.fig
%      IFEXTA, by itself, creates a new IFEXTA or raises the existing
%      singleton*.
%
%      H = IFEXTA returns the handle to a new IFEXTA or the handle to
%      the existing singleton*.
%
%      IFEXTA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IFEXTA.M with the given input arguments.
%
%      IFEXTA('Property','Value',...) creates a new IFEXTA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iFextA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iFextA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iFextA

% Last Modified by GUIDE v2.5 26-Sep-2017 15:38:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iFextA_OpeningFcn, ...
                   'gui_OutputFcn',  @iFextA_OutputFcn, ...
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


% --- Executes just before iFextA is made visible.
function iFextA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iFextA (see VARARGIN)

% Choose default command line output for iFextA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iFextA wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = iFextA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Cal_dsDNA.
function Cal_dsDNA_Callback(hObject, eventdata, handles)
% hObject    handle to Cal_dsDNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%dsDNAstruct=struct([]);
dsDNA_Px=[];
bp=str2num(get(handles.BP, 'string')); 
[dsr_filenm,dsr_path,dsr_Ind]=uigetfile('*.txt','Chose dsDNA Raw Data File');
RawDnaDataFile = strcat(dsr_path,dsr_filenm);
RawStgDat = importdata(RawDnaDataFile);

% [dsa_filenm,dsa_path,dsa_Ind] = uigetfile('*.txt','Chose dsDNA Analysed Data File');
% AnlDnaDataFile=strcat(dsa_path,dsa_filenm);
[~,nm,ext]=fileparts(dsr_filenm)
AnlDnaDataFile=sprintf('%s-analyzed.txt',nm);
AnlDnaDataFile= strcat(dsr_path,AnlDnaDataFile)

AnlStgDat = importdata(AnlDnaDataFile);
ImgIndx = RawStgDat(:,9);

% nonZIndx=find(ImgIndx);
% ImgIndx2=sort(ImgIndx(nonZIndx));
ind1=[];ind2=[];ind4=[];ind3=[];ind5=[];
datpath=uigetdir('','Chose dsDNA image directory');
dsDNA_Px = TrackExtPx( handles.Bead_1,handles.Bead_2,0.8,datpath,handles.Tracker,handles.I_R)

dsDNAstruct= struct('ImgInd1',ImgIndx,'rawStgExt',RawStgDat(:,2),'ImgInd2',dsDNA_Px(:,1),'ImgBead_2',dsDNA_Px(:,3),'ImgExt',dsDNA_Px(:,4),'AnalStr', AnlStgDat(:,2),'AnalStrImInd', AnlStgDat(:,4));

ind1=find(ismember(dsDNAstruct.ImgInd1,dsDNAstruct.ImgInd2));

nmExt=dsDNAstruct.rawStgExt(ind1);
pxExt=dsDNAstruct.ImgBead_2;
ind2=find(dsDNAstruct.ImgInd2<20000);
ind4=find(dsDNAstruct.ImgInd2>19999 & dsDNAstruct.ImgInd2<30000);
px2nmCal=mean(diff(nmExt(ind2)))/mean(abs(diff(pxExt(ind2))));

ind1=find(ismember(dsDNAstruct.ImgInd1,dsDNAstruct.ImgInd2));
ind3=find(ismember(dsDNAstruct.AnalStrImInd,dsDNAstruct.ImgInd2));
ind5=find(ismember(AnlStgDat(:,8),dsDNAstruct.ImgInd2));

dsExtstg=dsDNAstruct.AnalStr(ind3).*bp./1000;
dxExtIm=dsDNAstruct.ImgExt(ind2);

fun1= @(p) trapz(((px2nmCal.*dxExtIm-p)-dsExtstg).^2);
offset=fminsearch(fun1,3);

set(handles.Cal_Text,'string',num2str(px2nmCal));
set(handles.Offset_Text,'string',num2str(offset));

stgF_str=AnlStgDat(:,3);
stgF_rtn=AnlStgDat(:,7);
stgE_rtn=AnlStgDat(:,6);
imE_str=px2nmCal.*dsDNAstruct.ImgExt(ind2)-offset;
imE_str=(1000/bp).*imE_str;
imE_rtn=flipud(px2nmCal.*dsDNAstruct.ImgExt(ind4)-offset);
imE_rtn=(1000/bp).*imE_rtn;


axes(handles.MainAxes)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
h(1)=plot(dsDNAstruct.AnalStr,stgF_str,'b.','MarkerSize',15); hold on
h(3)=plot(stgE_rtn,stgF_rtn,'bo','MarkerSize',5)

h(2)=plot(imE_str,stgF_str(ind3),'r.','MarkerSize',15); hold on
h(4)=plot(imE_rtn,stgF_rtn(ind5),'ro','MarkerSize',5)
xlabel('Extension (nm/bp)','Color',[0.95 0.95 0.95])
ylabel('Force (pN)','Color',[0.95 0.95 0.95])
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
legend(h(1:2),'stage','image','Location','northwest')
handles.offset=offset;
handles.calib=px2nmCal;
guidata(hObject, handles);

% --- Executes on button press in Chose_ref_Beads.
function Chose_ref_Beads_Callback(hObject, eventdata, handles)
% hObject    handle to Chose_ref_Beads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Bead_1, Bead_2,I_R,I_ref]=Chose_Ref_Beads( );
axes(handles.Bead1)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5]);
imshow(Bead_1,[]);
axes(handles.Bead2)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5])
imshow(Bead_2,[]);

handles.Bead_1=Bead_1;
handles.Bead_2=Bead_2;
handles.I_R=I_R;

axes(handles.Tracker)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5])
imshow(I_ref,[]);
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TrackData.
function TrackData_Callback(hObject, eventdata, handles)
% hObject    handle to TrackData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str_im=[];rtn_im=[];ffe_im=[];ffr_im=[];str_stg=[];rtn_stg=[];ffe_stg=[];ffr_stg=[]; data_Px =[];ind1=[];ind2=[];ind3=[];ind4=[];ind5=[];ind6=[];
str_im_ind=[];
rtn_im_ind=[];
ffe_im_ind=[];
ffr_im_ind=[];
DataStruct= struct('Str_im_ind',str_im_ind,'Str_im',str_im,'Str_stg',str_stg,'Rtn_im_ind',rtn_im_ind,'Rtn_im',rtn_im,'Rtn_stg',rtn_stg,'Ffe_im_ind',ffe_im_ind,'Ffe_im',ffe_im,'Ffe_stg',ffe_stg,'Ffr_im_ind',ffr_im_ind','Ffr_im',ffr_im,'Ffr_stg',ffr_stg);

datpath=uigetdir('','Chose data image directory');
data_Px = TrackExtPx( handles.Bead_1,handles.Bead_2,0.8,datpath,handles.Tracker,handles.I_R);

data_Px(:,4)=handles.calib.*data_Px(:,4) - handles.offset; % Convert Pxl to um

ind1=find(data_Px(:,1)>9999 & data_Px(:,1)<20000); %find stretch images
ind2=find(data_Px(:,1)>19999 & data_Px(:,1)<30000); %find Return images
ind3=find(data_Px(:,1)>29999 & data_Px(:,1)<40000); %find Force feedback during extension images
ind4=find(data_Px(:,1)>39999);                      %find Force feedback during return images
[~,ind5]=sort(data_Px(ind2,4));


DataStruct.Str_im_ind= data_Px(ind1,1); DataStruct.Str_im= data_Px(ind1,4);
DataStruct.Rtn_im_ind= data_Px(ind2,1); DataStruct.Rtn_im= data_Px(ind2,4);
DataStruct.Rtn_im_ind=DataStruct.Rtn_im_ind(ind5) ; DataStruct.Rtn_im= DataStruct.Rtn_im(ind5);
DataStruct.Ffe_im_ind= data_Px(ind3,1); DataStruct.Ffe_im= data_Px(ind3,4);
DataStruct.Ffr_im_ind= data_Px(ind4,1); DataStruct.Ffr_im= data_Px(ind4,4);

handles.DataStructs=DataStruct;

% axes(handles.MainAxes)
% cla(handles.MainAxes)
% set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
% hold on
% [~,ind6]=sort(data_Px(:,1));
% plot(data_Px(ind6,1),data_Px(ind6,4),'b.','MarkerSize',15);
% xlabel('Image No.')
% ylabel('Extension (nm/bp)')
 guidata(hObject, handles);



function Cal_Text_Callback(hObject, eventdata, handles)
% hObject    handle to Cal_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cal_Text as text
%        str2double(get(hObject,'String')) returns contents of Cal_Text as a double


% --- Executes during object creation, after setting all properties.
function Cal_Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cal_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Offset_Text_Callback(hObject, eventdata, handles)
% hObject    handle to Offset_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Offset_Text as text
%        str2double(get(hObject,'String')) returns contents of Offset_Text as a double


% --- Executes during object creation, after setting all properties.
function Offset_Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Offset_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Import_stg_str.
function Import_stg_str_Callback(hObject, eventdata, handles)
% hObject    handle to Import_stg_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bp=str2num(get(handles.BP, 'string')); 
drift_str=[];
drift_rtn=[];
DataStruct=handles.DataStructs;
ind1=[];ind2=[];ind3=[];ind3=[];
[dsa_filenm,dsa_path,dsa_Ind] = uigetfile('*.txt','Chose Analyzed(stage) Stretch&Return Data File');
AnlDataFile=strcat(dsa_path,dsa_filenm);
AnlStgDat = importdata(AnlDataFile);
rtn_imInd=AnlStgDat(:,8);
str_imInd=AnlStgDat(:,4); 
ind1=find(ismember(DataStruct.Str_im_ind,str_imInd));
ind2=find(ismember(DataStruct.Rtn_im_ind,rtn_imInd));
t_str_stg=AnlStgDat(:,1); 
t_rtn_stg=AnlStgDat(:,5);
t_str_im=t_str_stg(ind1);
t_rtn_im=t_rtn_stg(ind2);
str_stg=AnlStgDat(:,2); 
str_im=DataStruct.Str_im.*(1000/bp);
str_F_stg=AnlStgDat(:,3);
str_F_im=str_F_stg(ind1);
rtn_stg=AnlStgDat(:,6); 
rtn_im=DataStruct.Rtn_im.*(1000/bp);
rtn_F_stg=AnlStgDat(:,7);
rtn_F_im=rtn_F_stg(ind2);
drift_str=str_stg(ind1)-str_im;
drift_rtn=rtn_stg(ind2)-rtn_im;






axes(handles.MainAxes)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
hold on
h(1)=plot(str_stg,str_F_stg,'b.','MarkerSize',15)
h(3)=plot(rtn_stg,rtn_F_stg,'bo','MarkerSize',5)
h(2)=plot(str_im,str_F_im,'r.','MarkerSize',15)
h(4)=plot(rtn_im,rtn_F_im,'ro','MarkerSize',5)
xlabel('Extension (nm/bp)','Color',[0.95 0.95 0.95])
ylabel('Force (pN)','Color',[0.95 0.95 0.95])
legend(h(1:2),'Stage','Image')


axes(handles.DriftAxes)
cla(handles.DriftAxes)
hold on
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
plot(t_str_im,drift_str,'g.','MarkerSize',8);
plot(t_rtn_im,drift_rtn,'go','MarkerSize',4);
xlabel('time (s)','Color',[0.95 0.95 0.95]);
ylabel('Extension (nm/bp)','Color',[0.95 0.95 0.95]);

handles.str_stg=str_stg; 
handles.str_im=str_im;
handles.str_F_stg=str_F_stg;
handles.str_F_im=str_F_im;
handles.rtn_stg=rtn_stg; 
handles.rtn_im=rtn_im;
handles.rtn_F_stg=rtn_F_stg;
handles.rtn_F_im=rtn_F_im;
handles.drift_str=drift_str;
handles.drift_rtn=drift_rtn;

handles.t_str_stg=t_str_stg; 
handles.t_rtn_stg=t_rtn_stg;
handles.t_str_im=t_str_im;
handles.t_rtn_im=t_rtn_im;


handles.AnlStrRtnStgData=AnlStgDat;
handles.datpath=dsa_path;
handles.datfilenm=dsa_filenm;

guidata(hObject, handles);
% DataStruct.Str_im_ind= data_Px(ind1,1); DataStruct.Str_im= data_Px(ind1,4);
% DataStruct.Rtn_im_ind= data_Px(ind2,1); DataStruct.Rtn_im= data_Px(ind2,4);
% DataStruct.Ffe_im_ind= data_Px(ind3,1); DataStruct.Ffe_im= data_Px(ind3,4);
% DataStruct.Ffr_im_ind= data_Px(ind4,1); DataStruct.Ffr_im= data_Px(ind4,4);


% --- Executes on button press in CorrectStrRtn.
function CorrectStrRtn_Callback(hObject, eventdata, handles)
DataStruct=handles.DataStructs;
drift_str_intrp =[];drift_str_intrp=[];
t_str_stg=handles.t_str_stg; 
t_rtn_stg=handles.t_rtn_stg;
t_str_im=handles.t_str_im;
t_rtn_im=handles.t_rtn_im;

str_stg=handles.str_stg; 
str_im=handles.str_im;
str_F_stg=handles.str_F_stg;
str_F_im=handles.str_F_im;
rtn_stg=handles.rtn_stg; 
rtn_im=handles.rtn_im;
rtn_F_stg=handles.rtn_F_stg;
rtn_F_im=handles.rtn_F_im;
drift_str=handles.drift_str;
drift_rtn=handles.drift_rtn;
AnlStgDat=handles.AnlStrRtnStgData;
dsa_path=handles.datpath;
dsa_filenm=handles.datfilenm;

MovAvgWin=sqrt(numel(drift_str));

drift_str_smth=movmean(drift_str,MovAvgWin);
drift_rtn_smth=movmean(drift_rtn,MovAvgWin)


drift_str_intrp = interp1(t_str_im,drift_str_smth,t_str_stg,'linear','extrap');
drift_rtn_intrp = interp1(t_rtn_im,drift_rtn_smth,t_rtn_stg,'linear','extrap');

str_stg_crctd=str_stg-drift_str_intrp;
rtn_stg_crctd=rtn_stg-drift_rtn_intrp;

axes(handles.MainAxes)
hold on
h(5)=plot(str_stg_crctd,str_F_stg,'m.','MarkerSize',15);
plot(rtn_stg_crctd,rtn_F_stg,'mo','MarkerSize',5)
legend(h(5),'Corrected Stage','Location','northwest')


axes(handles.DriftAxes)
cla(handles.DriftAxes)
hold on
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
plot(t_str_im,drift_str,'g.','MarkerSize',8);
plot(t_rtn_im,drift_rtn,'go','MarkerSize',4);
plot(t_str_stg,drift_str_intrp,'k-','MarkerSize',7);
plot(t_rtn_stg,drift_rtn_intrp,'k-','MarkerSize',7);
xlabel('time (s)','Color',[0.95 0.95 0.95]);
ylabel('Extension (nm/bp)','Color',[0.95 0.95 0.95]);



AnlStgDat=handles.AnlStrRtnStgData;
AnlStgDat(:,2)=str_stg_crctd;
AnlStgDat(:,6)=rtn_stg_crctd;
datpth=handles.datpath;
datfilenm=handles.datfilenm;
[~,nm,ext]=fileparts(datfilenm)
datfilenm=sprintf('%s_imCrctd.txt',nm);
datfilenm= strcat(datpth,datfilenm);

dlmwrite(datfilenm,AnlStgDat,'delimiter','\t','newline','pc');
% hObject    handle to CorrectStrRtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in imprt_FFE.
function imprt_FFE_Callback(hObject, eventdata, handles)
% hObject    handle to imprt_FFE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bp=str2num(get(handles.BP, 'string')); 
DataStruct=handles.DataStructs;
ind1=[];ind2=[];ind3=[];ind3=[];
[dsa_filenm,dsa_path,dsa_Ind] = uigetfile('*.txt','Chose Analyzed(stage) FFE Data File');
AnlDataFile=strcat(dsa_path,dsa_filenm);
AnlStgDat = importdata(AnlDataFile);


ffe_imgInd=AnlStgDat(:,4); 
ind1=find(ismember(ffe_imgInd,DataStruct.Ffe_im_ind));

t_stg=AnlStgDat(:,1); 
t_im=t_stg(ind1);
ext_stg=AnlStgDat(:,2);
F_stg=AnlStgDat(:,3); 
ext_im=DataStruct.Ffe_im.*(1000/bp);
F_im=F_stg(ind1);
drift=ext_stg(ind1)-ext_im;

handles.ffe_ext_im=ext_im;
handles.ffe_ext_stg=ext_stg;
handles.ffe_t_stg=t_stg;
handles.ffe_t_im=t_im;
handles.AnlFFEData=AnlStgDat;
handles.datFFEpath=dsa_path;
handles.datFFEfilenm=dsa_filenm;


axes(handles.MainAxes)
cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
hold on
h(1)=plot(t_stg,ext_stg,'b.','MarkerSize',15);
h(2)=plot(t_im,ext_im,'r.','MarkerSize',15);

xlabel('time (s)','Color',[0.95 0.95 0.95]);
ylabel('Extension (nm/bp)','Color',[0.95 0.95 0.95]);
legend(h(1:2),'Stage','Image')

axes(handles.DriftAxes)
cla(handles.DriftAxes)
hold on
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
plot(t_im,drift,'y.','MarkerSize',7);
xlabel('time (s)','Color',[0.95 0.95 0.95]);
ylabel('Extension (nm/bp)','Color',[0.95 0.95 0.95]);

handles.ffe_drift=drift;
handles.ffe_ext_im=ext_im;
handles.ffe_ext_stg=ext_stg;
handles.ffe_t_stg=t_stg;
handles.ffe_t_im=t_im;
handles.AnlFFEData=AnlStgDat;
handles.datFFEpath=dsa_path;
handles.datFFEfilenm=dsa_filenm;

guidata(hObject, handles);

% --- Executes on button press in CrectFFE.
function CrectFFE_Callback(hObject, eventdata, handles)
% hObject    handle to CrectFFE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ext_im=handles.ffe_ext_im;
ext_stg=handles.ffe_ext_stg;
t_stg=handles.ffe_t_stg;
t_im=handles.ffe_t_im;
AnlStgDat=handles.AnlFFEData;
datpth=handles.datFFEpath;
datfilenm=handles.datFFEfilenm;
drift=handles.ffe_drift;

MovAvgWin=sqrt(numel(drift));
drift_Smth=movmean(drift,MovAvgWin);
drift_interp = interp1(t_im,drift_Smth,t_stg,'linear','extrap');
ext_stg_crctd=ext_stg-drift_interp;

cla(handles.DriftAxes)
hold on
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
plot(t_im,drift,'y.','MarkerSize',7);
plot(t_im,drift_Smth,'k-','LineWidth',2);
xlabel('time (s)','Color',[0.95 0.95 0.95]);
ylabel('Extension (nm/bp)','Color',[0.95 0.95 0.95]);

axes(handles.MainAxes)
%cla(handles.MainAxes)
set(gca,'Color',[0.5 0.5 0.5],'XColor',[0.95 0.95 0.95],'YColor',[0.95 0.95 0.95]);
hold on
h(1)=plot(t_stg,ext_stg_crctd,'g.','MarkerSize',15);
legend(h(1),'Corrected Stage data');

AnlStgDat(:,2)=ext_stg_crctd;

[~,nm,ext]=fileparts(datfilenm)
datfilenm=sprintf('%s_imCrctd.txt',nm);
datfilenm= strcat(datpth,datfilenm);
dlmwrite(datfilenm,AnlStgDat,'delimiter','\t','newline','pc');

% --- Executes on button press in ImprtFFR.
function ImprtFFR_Callback(hObject, eventdata, handles)
% hObject    handle to ImprtFFR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CrectFFR.
function CrectFFR_Callback(hObject, eventdata, handles)
% hObject    handle to CrectFFR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BP_Callback(hObject, eventdata, handles)
% hObject    handle to BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BP as text
%        str2double(get(hObject,'String')) returns contents of BP as a double


% --- Executes during object creation, after setting all properties.
function BP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

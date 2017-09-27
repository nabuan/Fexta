function varargout = gui(varargin)


% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 27-Oct-2016 15:29:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_OpeningFcn, ...
    'gui_OutputFcn',  @gui_OutputFcn, ...
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

% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

stiffflag=0;
handles.stiffflag=stiffflag;
handles.DataFileflag=0;
guidata(hObject, handles);
handles.FilterIndexData=0;
handles.DNAanalysedFlag=0;
guidata(hObject, handles);
cla(handles.DNA)


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function pushbutton6_Callback(hObject, eventdata, handles)
handles.orgdir=pwd;
FilterIndexStiff=0;
[Fst,stpath,FilterIndexStiff] = uigetfile({'*.*','All Files(*.*)'}, 'provide the stiffness file ');

if(FilterIndexStiff==1)
    st = strcat(stpath,Fst);
    stiffFile=importdata(st);
    
    handles.stiffflag;
    stiffflag=1;
    handles.stiffflag=stiffflag;
    handles.stiffFile=stiffFile;
    handles.newdir=stpath;
    set(handles.stifftick,'value',1);
    set(handles.display,'string','stiffness file was loaded');
    handles.DNAanalysedFlag=0;
    guidata(hObject, handles);
else
    
    warndlg('Failed to load stiffness file')
end


% --- Executes on button press in pushbutton6.
function AnaDNA_Callback(hObject, eventdata, handles)
% % --- Executes on button press in AnaDNA.
stiffflag=handles.stiffflag;
if(stiffflag ==0)
    %
    warndlg('Please Load stiffness file')
else
    
    cd(handles.newdir)
    
    [F,PathName,fileFlagDNA] = uigetfile({'*.*','All Files(*.*)'}, 'Select your DNA File ');
    
    handles.newdir=PathName;
    if(fileFlagDNA==1)
        handles.FilterIndexData=0;
        cla(handles.DNA)
        
        b = strcat(PathName,F);
        
        DNA=importdata(b);
        n=find(all(DNA==0,2));
        for(i=1:length(n))
            jj=n(i);
            if(all(DNA(jj+1,:)==0))
                DNA(jj+1,:)=[];
                break;
            end
        end
        
        
        n=n(1);
        teDNA(1,:)=DNA((1:n-1),1);
        xeDNA(1,:)=DNA((1:n-1),2);
        feDNA(1,:)=DNA((1:n-1),3);
        feDNA(2,:)=DNA((1:n-1),6);
        
        trDNA(1,:)=DNA((n+1:end),1);
        xrDNA(1,:)=DNA((n+1:end),2);
        frDNA(1,:)=DNA((n+1:end),3);
        frDNA(2,:)=DNA((n+1:end),6);
        
        %write a condition to check if this column exist
        img_ind_e(1,:)=DNA((1:n-1),9);
        img_ind_r(1,:)=DNA((n+1:end),9);
        
        handles.ret=xrDNA;
        handles.Fret=frDNA;
        
        feDNARaw=mean(feDNA);
        Fexten=feDNA;
        exten=xeDNA;
        handles.feDNARaw=feDNARaw;
        handles.teDNA=teDNA;
        handles.RawXeDNA=xeDNA;
        handles.RawXrDNA=xrDNA;
        handles.RawFeDNA=feDNA;
        handles.RawFrDNA=frDNA;
        handles.DataFileflag=0;
        handles.PathName=PathName;
        handles.FilePartsDNA=F;
        handles.ext = xeDNA;
        handles.Fext1=feDNA(1,:);
        handles.Fext2=feDNA(2,:);
        guidata(hObject, handles);
        
        tol_str=get(handles.tolerance,'string');
        tol=str2double(tol_str);
        stiffOrder=get(handles.stiffOrderstr,'val')
        % stiffOrder=str2double(stiffOrder_str)
        bp_str=get(handles.bp, 'string');
        Salt_str=get(handles.midpoint, 'string');
        SaltConc=str2double(Salt_str);
        if(SaltConc>250 || SaltConc<2.6)
            warndlg('Salt concentration is out of acceptable range')
            cd(handles.orgdir);
            return;
        else
            cd(handles.orgdir);
            Fmp=ForceMidpoint(SaltConc);
            cd(handles.PathName);
            FmpString=sprintf('Force Midpoint at %0.2f pN',Fmp);
            set(handles.display,'string',FmpString);
        end
        bp= str2double(bp_str);
        stiffFile = handles.stiffFile;
        cd(handles.orgdir);
        Fmp=ForceMidpoint(SaltConc);
        
        [calibration,xoffset,blWLC,bL,stiffCoeff,EWLC,FWLC,minplat,maxplat]=main(bp,Fmp,feDNA,xeDNA,stiffFile,tol,stiffOrder);
        [FeDNA,xeDNA]=AnalayzeData(bp,feDNA,xeDNA,calibration,blWLC,bL,stiffCoeff,xoffset);
        
        set(handles.BL01,'string',num2str(bL(1)));
        set(handles.BL02,'string',num2str(bL(2)));
        set(handles.bLwlc,'string',num2str(blWLC));
        set(handles.xOffset,'string',num2str(xoffset));
        set(handles.cAlibration,'string',num2str(prod(calibration)));
        set(handles.sTiffness,'string',num2str(1/stiffCoeff(stiffOrder)));
        
        axes(handles.DNA);
        
        plot(xeDNA,FeDNA,'ob','MarkerSize',4,'MarkerFaceColor','b');hold on
        line([xeDNA(minplat),xeDNA(minplat)],[min(FeDNA),max(FeDNA)],'color','g','LineStyle',':','linewidth',3)
        line([xeDNA(maxplat),xeDNA(maxplat)],[min(FeDNA),max(FeDNA)],'color','g','LineStyle',':','linewidth',3)
        line([min(xeDNA),max(xeDNA)],[Fmp,Fmp],'color','m','LineStyle',':','linewidth',3)
        xlabel('Extension (nm/bp)')
        ylabel('Force (pN)')
        
        grid on;
        
        [FrDNA,xrDNA]=AnalayzeData(bp,frDNA,xrDNA,calibration,blWLC,bL,stiffCoeff,xoffset);
        handles.calibration=calibration;
        handles.calibrationORG=calibration;
        handles.stiffOrderORG=stiffOrder;
        handles.xoffset=xoffset;
        handles.xoffsetORG=xoffset;
        handles.baseline=blWLC;
        handles.baselineORG=blWLC;
        handles.bL=bL;
        handles.stiffCoeff=stiffCoeff;
        handles.stiffCoeffORG=stiffCoeff;
        handles.DNAanalysedFlag=1;
        % handles.stiffness=stiffness;
        % handles.stiffnessORG=stiffness;
        
        %%%%%%%pulling rate%%%%%%%%%%%%%%%%
        pullingRate=mean( diff(bp.*xeDNA)./diff(teDNA));
        pullingRateErr=std(diff(bp.*xeDNA)./diff(teDNA))/sqrt(length(diff(teDNA)));
        pullRateString=sprintf('%0.2f +/- %0.2f',pullingRate,pullingRateErr);
        set(handles.pullRate,'string',pullRateString);
        
        [xrDNA,i]=sort(xrDNA);
        FrDNA=FrDNA(i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.xeDNA=xeDNA;
        handles.xrDNA=xrDNA;
        handles.FeDNA=FeDNA;
        handles.FrDNA=FrDNA;
        handles.basep=bp;
        %handles.StiffnessData=StiffnessData;
        handles.EWLC=EWLC;
        handles.FWLC=FWLC;
        handles.minplat=minplat;
        handles.maxplat=maxplat;
        handles.Fmp=Fmp;
        guidata(hObject, handles);
        
        
        axes(handles.DNA);
        plot(xrDNA,FrDNA,'ob','MarkerSize',4);
        plot(EWLC,FWLC,'r','linewidth',2)
        
        ymax=max(max(handles.FeDNA))+5;
        ymin=0;
        xmax=max(xeDNA)+0.02;
        xmin=min(xeDNA)-0.01;
        axis([xmin xmax ymin ymax]);
        set(handles.DNAtick,'value',1);
        
        cd(handles.orgdir)
    else
        warndlg('Failed to Load DNA only file')
        cd(handles.orgdir)
        
    end
    %%%%%%%%%%%%%Auto Saving%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % xeDNA=handles.xeDNA;
    % xrDNA=handles.xrDNA;
    % FeDNA=handles.FeDNA;
    % FrDNA=handles.FrDNA;
    
    
    AutoSaveVal = get(handles.AutoSave, 'Value')
    if(AutoSaveVal==1)
        DNAonly=handles.FilePartsDNA;
        [pathstr,DNAname,ext] = fileparts(DNAonly);
        analyzedDNA= sprintf('%s-analyzed.txt',DNAname);
        M=reshapeData([teDNA;xeDNA;FeDNA;img_ind_e],[trDNA;xrDNA;FrDNA;img_ind_r]);
        cd(handles.PathName);
        dlmwrite(analyzedDNA,M,'delimiter','\t','newline','pc')
        set(handles.display,'string','Analyzed DNA Only File was saved');
        cd(handles.orgdir)
        
    end
end


function Legend_Callback(hObject, eventdata, handles)
% DNA only
% Data Only
% DNA and Data
% Stiffness
% Raw Data - DNA only
% Raw Data - Data
% Raw Data - Stiffness
if(handles.DNAanalysedFlag==1)
    xeDNA=handles.xeDNA;
    xrDNA=handles.xrDNA;
    FeDNA=handles.FeDNA;
    FrDNA=handles.FrDNA;
end
%StiffnessData=handles.StiffnessData;

axes(handles.DNA);
val = get(hObject,'Value');


switch val
    case 1
        cla(handles.DNA)
        display('DNAonly');
        plot(xeDNA,FeDNA,'ob','MarkerSize',4,'MarkerFaceColor','b');hold on
        plot(xrDNA,FrDNA,'ob','MarkerSize',4)
        
        ymax=max(handles.FeDNA)+5;
        ymin=0;
        xmax=max(xeDNA)+0.02;
        xmax=max(xeDNA)+0.02;
        xmin=min(xeDNA)-0.01;
        axis([xmin xmax ymin ymax]);
        xlabel('Extension (nm/bp)')
        ylabel('Force (pN)')
    case 3
        FeData=handles.FeData;
        FrData=handles.FrData;
        xeData=handles.xeData;
        xrData=handles.xrData;
        
        cla(handles.DNA)
        display('DNA and Data');
        plot(xeDNA,FeDNA,'ob','MarkerSize',4,'MarkerFaceColor','b');hold on
        plot(xrDNA,FrDNA,'ob','MarkerSize',4)
        plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
        plot(xrData,FrData,'or','MarkerSize',4)
        ymax=max(handles.FeDNA)+5;
        ymin=0;
        xmax=max(xeDNA)+0.02;
        xmax=max(xeDNA)+0.02;
        xmin=min(xeDNA)-0.01;
        axis([xmin xmax ymin ymax]);
        xlabel('Extension (nm/bp)')
        ylabel('Force (pN)')
    case 2
        FeData=handles.FeData;
        FrData=handles.FrData;
        xeData=handles.xeData;
        xrData=handles.xrData;
        
        cla(handles.DNA)
        display('DNA and Data');
        plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
        plot(xrData,FrData,'or','MarkerSize',4)
        ymax=max(handles.FeData)+5;
        ymin=0;
        xmax=max(xeData)+0.02;
        xmax=max(xeData)+0.02;
        xmin=min(xeData)-0.01;
        axis([xmin xmax ymin ymax]);
        xlabel('Extension (nm/bp)')
        ylabel('Force (pN)')
    case 4
        
        
        calibration=handles.calibration;
        
        
        stiffFile=handles.stiffFile;
        cla(handles.DNA)
        n=find(all(stiffFile==0,2));
        Stext(1,:)=stiffFile((1:n-1),2);
        StFext(1,:)=stiffFile((1:n-1),3);
        StFext(2,:)=stiffFile((1:n-1),6);
        StFext=calibration.*mean(StFext);
        stiffOrder=get(handles.stiffOrderstr,'val');
        stiffCoeff= findStiffnessFit(calibration,stiffFile,stiffOrder);
        XX=correctStiffness(stiffCoeff,Stext,StFext,1000,1);
        plot(XX,StFext,'ob','MarkerSize',4,'MarkerFaceColor','b');hold on
        ymax=max(StFext)+ 0.1*max(StFext);
        ymin=min(StFext)- 0.1*max(StFext);
        xmax=max(XX)+1;
        xmin=min(XX)-1;
        xlabel('extension(um)')
        ylabel('Force (pN)')
        axis([xmin xmax ymin ymax]);
        handles.stiffCoeff=stiffCoeff;
        guidata(hObject,handles);
        
    case 5
        
        % xrDNA=handles.RawXrDNA;
        % Fret1=handles.RawFeDNA(1,:);
        % Fret2=handles.RawFeDNA(2,:);
        % [xrDNA,i]=sort(xrDNA)
        % Fret1=Fret1(i)
        % Fret2=Fret2(i)
        
        cla(handles.DNA)
        plot(handles.ext,handles.Fext1,'ob','MarkerSize',4,'MarkerFaceColor','b'); hold on
        plot(handles.ext,handles.Fext2,'og','MarkerSize',4,'MarkerFaceColor','g');
        plot(handles.ret,handles.Fret(1,:),'ob','MarkerSize',4); hold on
        plot(handles.ret,handles.Fret(2,:),'og','MarkerSize',4);
        ymax=max(max(handles.Fext1,handles.Fext2))+0.1*max(max(handles.Fext1,handles.Fext2));
        ymin=min(min(handles.Fext1,handles.Fext2));
        xmax=max(max(handles.ext))+1;
        axis([0 xmax ymin ymax]);
        xlabel('extension(um)')
        ylabel('Force (au)')
    case 6
        if(handles.DataFileflag==1);
            Data=handles.RawData ;
            tData(1,:)=Data(:,1);
            xData(1,:)=Data(:,2);
            FData(1,:)=Data(:,3);
            FData(2,:)=Data(:,6);
            cla(handles.DNA);
            plot(xData,FData(1,:),'ob','MarkerSize',4,'MarkerFaceColor','b'); hold on
            plot(xData,FData(2,:),'og','MarkerSize',4,'MarkerFaceColor','g');
            ymax=max(max(FData(1,:),FData(2,:)))+0.1*max(max(FData(1,:),FData(2,:)));
            
            ymin=min(min(FData(1,:),FData(2,:)));
            xmax=max(max(xData))+1;
            axis([0 xmax ymin ymax]);
            xlabel('extension(um)')
            ylabel('Force (au)')
        end
    case 7
        if(handles.stiffflag==1)
            stiffFile=handles.stiffFile;
            n=find(all(stiffFile==0,2));
            xDataSt(1,:)=stiffFile((1:n-1),2);
            FextSt(1,:)=stiffFile((1:n-1),3);
            FextSt(2,:)=stiffFile((1:n-1),6);
            cla(handles.DNA)
            plot(xDataSt,FextSt(1,:),'ob','MarkerSize',4,'MarkerFaceColor','b'); hold on
            plot(xDataSt(1,:),FextSt(2,:),'og','MarkerSize',4,'MarkerFaceColor','g');
            ymax=max(max(FextSt))+10;
            ymin=min(min(FextSt))-10;
            xmax=max(max(xDataSt))+1;
            axis([0 xmax ymin ymax]);
            xlabel('extension(um)')
            ylabel('Force (au)')
        end
        
end




% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)

cla(handles.DNA)
cla(handles.FFgraph)
clc;
handles.stiffflag=0;
handles.DNAanalysedFlag=0;
set(handles.stifftick,'value',0);
set(handles.DNAtick,'value',0);
handles.FilterIndexData=0;
set(handles.display,'string','');
set(handles.display2,'string','');
set(handles.BL01,'string',num2str(0));
set(handles.sTiffness,'string',num2str(0));
set(handles.BL02,'string',num2str(0));
set(handles.bLwlc,'string',num2str(0));
set(handles.xOffset,'string',num2str(0));
set(handles.cAlibration,'string',num2str(prod(0)));

guidata(hObject, handles);

handles.stiffFile=[];
handles.RawData=[];
handles.stiffCoeff=[];
handles.DNAanalysedFlag=0;
H1=handles.H1;
H2=handles.H2;
delete(H2)
%delete(H1);
guidata(hObject, handles);
clear all


function AnalyzeDATA_Callback(hObject, eventdata, handles)
if(handles.DNAanalysedFlag==1 )
    cd(handles.newdir);
    [F2,PathNameData,FilterIndexData] = uigetfile({'*.*','All Files(*.*)'}, 'Select your DATA File ');
    handles.FilterIndexData=FilterIndexData;
    if(FilterIndexData==1)
        [pathstr,DataName,ext] = fileparts(F2);
        set(handles.DataFileName,'string',DataName);
        % --- Executes on button press in AnalyzeDATA.
        cla(handles.DNA)
        bp=handles.basep;
        % calibration=handles.calibration;
        % xoffset=handles.xoffset;
        % baseline=handles.baseline;
        % bL=handles.bL;
        % stiffness=handles.stiffness;
        xeDNA=handles.xeDNA;
        xrDNA=handles.xrDNA;
        FeDNA=handles.FeDNA;
        FrDNA=handles.FrDNA;
        
        calibration_str=get(handles.cAlibration,'string');
        calibration=str2double(calibration_str);
        xoffset_str=get(handles.xOffset,'string');
        xoffset=str2double(xoffset_str);
        baseline_str=get(handles.bLwlc,'string');
        baseline=str2double(baseline_str);
        bL01_str=get(handles.BL01,'string');
        bL(1)=str2double(bL01_str);
        bL02_str=get(handles.BL02,'string');
        bL(2)=str2double(bL02_str);
        stiffCoeff=handles.stiffCoeff;
        
        
        
        
        dataFile = strcat(PathNameData,F2);
        Data=importdata(dataFile);
        m=find(all(Data==0,2));
        
        for(i=1:length(m))
            jj=m(i);
            if(all(Data(jj+1,:)==0))
                Data(jj+1,:)=[];
                break;
            end
        end
        
        
        m=find(all(Data==0,2));
        type=length(m);
        k=length(Data);
        cd(handles.orgdir);
        
        tData(1,:)=Data(:,1);
        xData(1,:)=Data(:,2);
        FData(1,:)=Data(:,3);
        FData(2,:)=Data(:,6);
        
        imgIndData(1,:)=Data(:,9);
        
        handles.RawData=Data;
        
        cd(handles.orgdir);
        [FAData,xAData]=AnalayzeData(bp,FData,xData,calibration,baseline,bL,stiffCoeff,xoffset);
        
        handles.FAData=FAData;
        handles.xAData=xAData;
        handles.tData=tData;
        
        switch(type)
            case 1
                Extension=1:1:m(1)-1;
                Return=m(1)+1:1:k;
                set(handles.display2, 'String', 'detected only extension and return');
                teData=tData(Extension);
                xePdata=xAData(Extension); % to plot
            case 2
                Extension1=1:1:m(1)-1;
                teData=tData(Extension1);
                xePdata=xAData(Extension1);
                FF_Extension=m(1)+1:1:m(2)-1;
                Extension=[Extension1';FF_Extension'];
                Return=m(2)+1:1:k;
                
                cla(handles.FFgraph);
                axes(handles.FFgraph);
                [AX H1 H2]= plotyy(tData(FF_Extension),xAData(FF_Extension),tData(FF_Extension),FAData(FF_Extension));
                set(AX(2),'Ylim',[0 max(FAData(FF_Extension))+5],'YTick',[0:10:max(FAData(FF_Extension))+5],'YColor',[0.2,1,0.4],'XColor','w');
                set(AX(1),'Ylim',[min(xAData(FF_Extension))-0.01, max(xAData(FF_Extension))+0.01],'YTick',[0:0.01:max(xAData(FF_Extension))+0.01],'box','off','YColor',[0,0.8,1],'XColor','w');
                
                set(handles.display2, 'String', 'detected Forcefeedback on extension');
                ylabel(AX(2),'Force (pN)')
                ylabel(AX(1),'Extension (nm/bp)')
                xlabel('time (s)')
                handles.H1=H1;
                handles.H2=H2;
                handles.AX=AX;
                handles.AX1=AX(1);
                handles.AX2=AX(2);
                handles.FF_Extension=FF_Extension';
                handles.ExtensionNrm=Extension1';
                set(handles.FFmenu,'Value',1);
                guidata(hObject, handles);
                
            case 3
                Extension=1:1:m(1)-1;
                teData=tData(Extension);
                xePdata=xAData(Extension);
                Return1=m(1)+1:1:m(2)-1;
                FF_Return=m(2)+1:1:m(3)-1;
                Return2=m(2)+1:1:m(3)-1;
                Return3=m(3)+1:1:k;
                Return=[Return1';Return2';Return3'];
                axes(handles.FFgraph);
                cla(handles.FFgraph);
                [AX H1 H2]= plotyy(tData(FF_Return),xAData(FF_Return),tData(FF_Return),FAData(FF_Return));
                set(AX(2),'Ylim',[0 max(FAData(FF_Return))+5],'YTick',[0:10:max(FAData(FF_Return))+5],'YColor',[0.2,1,0.4]);
                set(AX(1),'Ylim',[min(xAData(FF_Return))-0.01 max(xAData(FF_Return))+0.01],'YTick',[0:0.01:max(xAData(FF_Return))+0.01],'box','off','YColor',[0,0.8,1],'XColor','w');
                
                ylabel(AX(2),'Force (pN)')
                ylabel(AX(1),'Extension (nm/bp)')
                xlabel('time (s)')
                
                set(handles.display2, 'String', 'detected Forcefeedback on return');
                handles.H1=H1
                handles.H2=H2
                handles.AX=AX;
                handles.AX1=AX(1);
                handles.AX2=AX(2);
                handles.FF_Return=FF_Return';
                handles.ExtensionNrm=Extension';
                handles.ReturnNrm=[Return1';Return3'];
                set(handles.FFmenu,'Value',2);
                guidata(hObject, handles);
            case 4
                Extension1=1:1:m(1)-1;
                teData=tData(Extension1);
                FF_Extension=m(1)+1:1:m(2)-1;
                Extension=[Extension1';FF_Extension'];
                xePdata=xAData(Extension1);
                Return1=m(2)+1:1:m(3)-1;
                FF_Return=m(3)+1:1:m(4)-1;
                Return2=FF_Return;
                Return3=m(4)+1:1:k;
                Return=[Return1';Return2';Return3'];
                axes(handles.FFgraph);
                [AX H1 H2]= plotyy(tData(FF_Extension),xAData(FF_Extension),tData(FF_Extension),FAData(FF_Extension));
                set(AX(2),'Ylim',[0 max(FAData(FF_Return))+5],'YTick',[0:10:max(FAData(FF_Return))+5],'YColor',[0.2,1,0.4],'XColor','w');
                set(AX(1),'Ylim',[min(xAData(FF_Return))-0.01 max(xAData(FF_Return))+0.01],'YTick',[0:0.01:max(xAData(FF_Return))+0.01],'box','off','YColor',[0,0.8,1],'XColor','w');
                
                ylabel(AX(2),'Force (pN)')
                ylabel(AX(1),'Extension (nm/bp)')
                xlabel('time (s)')
                
                handles.H1=H1
                handles.H2=H2
                handles.AX=AX;
                handles.AX1=AX(1);
                handles.AX2=AX(2);
                set(AX(2),'Ylim',[0 70]);
                set(handles.display2, 'String', 'detected Forcefeedback on return and extension');
                handles.ExtensionNrm=Extension';
                handles.FF_Return=FF_Return';
                handles.ExtensionNrm=Extension';
                handles.ReturnNrm=[Return1';Return3'];
                set(handles.FFmenu,'Value',2);
                guidata(hObject, handles);
        end
        
        xeData=xAData(Extension);
        xrData=xAData(Return);
        FeData(1,:)=FAData(Extension);
        FrData(1,:)=FAData(Return);
        
        handles.Extension=Extension;
        handles.teData=teData;
        handles.Return=Return;
        handles.xrData=xrData;
        handles.FrData=FrData;
        handles.xeData=xeData;
        handles.FeData=FeData;
        handles.DataFileflag=1;
        handles.DataType=type;
        handles.ZeroRows=m;
        handles.PathNameData=PathNameData;
        handles.filepartsData=F2;
        guidata(hObject, handles);
        axes(handles.DNA);
        
        plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
        plot(xrData,FrData,'or','MarkerSize',4)
        xlabel('Extension (nm/bp)')
        ylabel('Force (pN)')
        
        %%%%%%%%%%pullingRate%%%%%%%%%%%%%%%
        pullingRate=mean( diff(bp.*xePdata)./diff(teData));
        pullingRateErr=std( diff(bp.*xePdata)./diff(teData))/sqrt(length(diff(teData)));
        pullRateString=sprintf('%0.2f +/- %0.2f',pullingRate,pullingRateErr);
        set(handles.pullRate,'string',pullRateString);
        %%%%%%%%%pullingRate%%%%%%%%%%%%%%%%
        
        [xrDNA,i]=sort(xrDNA);
        FrDNA=FrDNA(i);
        
        AutoSaveVal = get(handles.AutoSave, 'Value');
        if(AutoSaveVal==1)
            [pathstr,DataName,ext] = fileparts(F2);
            [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,imgIndData,handles.ZeroRows);
            StrRertnFileName=sprintf('%s-anlyzd.txt',DataName);
            cd(PathNameData);
            dlmwrite(StrRertnFileName,StrRetnData,'delimiter','\t','newline','pc');
            handles.StrRetnData=StrRetnData;
            handles.tData_FFE=tData_FFE;
            handles.tData_FFR=tData_FFR;
            guidata(hObject, handles);
            switch(handles.DataType)
                case 2
                    FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                    dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                case 3
                    FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                    dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
                case 4
                    FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                    dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                    FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                    dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
            end
            cd(handles.orgdir)
            set(handles.display,'string','Analyzed Data File(s) was/were saved');
        end
    else
        warndlg('Failed to Load Data file')
        cd(handles.orgdir)
        set(handles.Legend,'Value',2);
    end
    
else
    warndlg('Couldnt detect analyzed DNA')
end
% --------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function pushbutton6_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to pushbutton6 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
%
% % --- Executes on button press in pushbutton6.


% --- Executes on button press in gRid.
function gRid_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if(val==1)
    axes(handles.DNA);
    grid on
else
    grid off
end


% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
stiffFile=handles.stiffFile;
bp=handles.basep;
calibration_str=get(handles.cAlibration,'string');
calibration=str2double(calibration_str);
handles.calibration=calibration;
xoffset_str=get(handles.xOffset,'string');
xoffset=str2double(xoffset_str);
baseline_str=get(handles.bLwlc,'string');
baseline=str2double(baseline_str);
bL01_str=get(handles.BL01,'string');
bL(1)=str2double(bL01_str);
bL02_str=get(handles.BL02,'string');
bL(2)=str2double(bL02_str);
% stiffness_str=get(handles.sTiffness,'string');
% stiffness=str2double(stiffness_str);
stiffOrder=get(handles.stiffOrderstr,'val');
stiffCoeff= findStiffnessFit(calibration,stiffFile,stiffOrder);
set(handles.sTiffness,'string',num2str(1/stiffCoeff(stiffOrder)));
if(handles.FilterIndexData==1)
    Extension=handles.Extension;
    Return=handles.Return;
    
    
    Data=handles.RawData;
    
    tData(1,:)=Data(:,1);
    xData(1,:)=Data(:,2);
    FData(1,:)=Data(:,3);
    FData(2,:)=Data(:,6);
    
    
    
    [FAData,xAData]=AnalayzeData(bp,FData,xData,calibration,baseline,bL,stiffCoeff,xoffset);
    
    handles.FAData=FAData;
    handles.xAData=xAData;
    handles.tData=tData;
    xeData=xAData(Extension);
    xrData=xAData(Return);
    FeData(1,:)=FAData(Extension);
    FrData(1,:)=FAData(Return);
    
    handles.xrData=xrData;
    handles.FrData=FrData;
    handles.xeData=xeData;
    handles.FeData=FeData;
    
    % handles.stiffness=stiffness;
    handles.calibration=calibration;
    guidata(hObject, handles);
    axes(handles.DNA);
    cla(handles.DNA)
    plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
    plot(xrData,FrData,'or','MarkerSize',4)
    ymax=max(handles.FeData)+5;
    ymin=0;
    xmax=max(xeData)+0.02;
    xmax=max(xeData)+0.02;
    xmin=min(xeData)-0.01;
    axis([xmin xmax ymin ymax]);
    [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,handles.ZeroRows);
    handles.StrRetnData=StrRetnData;
    handles.tData_FFE=tData_FFE;
    handles.tData_FFR=tData_FFR;
    guidata(hObject, handles);
    AutoSaveVal = get(handles.AutoSave, 'Value');
    if(AutoSaveVal==1)
        [pathstr,DataName,ext] = fileparts(handles.filepartsData);
        StrRertnFileName=sprintf('%s-anlyzd.txt',DataName);
        cd(handles.PathNameData);
        dlmwrite(StrRertnFileName,StrRetnData,'delimiter','\t','newline','pc');
        set(handles.display,'string','Analyzed Data File was updated');
    end
    
    
    switch(handles.DataType)
        case 2
            axes(handles.FFgraph)
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFE(:,1),tData_FFE(:,2),tData_FFE(:,1),tData_FFE(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFE(:,3))+5],'YTick',[0:10:max(tData_FFE(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFE(:,2))-0.01 max(tData_FFE(:,2))+0.01],'YTick',[0:0.01:max(tData_FFE(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                set(handles.display,'string','Analyzed Data File was updated');
            end
        case 3
            
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            if(AutoSaveVal==1)
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
            end
        case 4
            axes(handles.FFgraph)
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
                set(handles.display,'string','Analyzed Data File was updated');
            end
    end
    cd(handles.orgdir)
else
    xeDNA=handles.RawXeDNA;
    xrDNA=handles.RawXrDNA;
    feDNA=handles.RawFeDNA;
    frDNA=handles.RawFrDNA;
    
    
    
    [FrDNA,xrDNA]=AnalayzeData(bp,frDNA,xrDNA,calibration,baseline,bL,stiffCoeff,xoffset);
    [FeDNA,xeDNA]=AnalayzeData(bp,feDNA,xeDNA,calibration,baseline,bL,stiffCoeff,xoffset);
    
    cla(handles.DNA)
    axes(handles.DNA)
    plot(xrDNA,FrDNA,'ob','MarkerSize',4);
    plot(xeDNA,FeDNA,'ob','MarkerSize',4,'MarkerFaceColor','b');
    
    ymax=max(max(handles.FeDNA))+5;
    ymin=0;
    xmax=max(xeDNA)+0.02;
    
    xmin=min(xeDNA)-0.01;
    axis([xmin xmax ymin ymax]);
    set(handles.DNAtick,'value',1);
    
    AutoSaveVal = get(handles.AutoSave, 'Value')
    [xrDNA i]=sort(xrDNA);
    FrDNA=FrDNA(i);
    
    handles.xeDNA=xeDNA;
    handles.xrDNA=xrDNA;
    handles.FeDNA=FeDNA;
    handles.FrDNA=FrDNA;
    guidata(hObject, handles);
    
    if(AutoSaveVal==1)
        DNAonly=handles.FilePartsDNA;
        [pathstr,DNAname,ext] = fileparts(handles.FilePartsDNA);
        analyzedDNA= sprintf('%s-analyzed.txt',DNAname);
        M=reshapeData([xeDNA;FeDNA],[xrDNA;FrDNA]);
        cd(handles.PathName);
        dlmwrite(analyzedDNA,M,'delimiter','\t','newline','pc')
        set(handles.display,'string','Analyzed DNA Only File was updated');
        cd(handles.orgdir)
    end
    cd(handles.orgdir)
end


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)

stiffOrder=handles.stiffOrderORG;

stiffCoeff=handles.stiffCoeffORG;
set(handles.sTiffness,'string',num2str(1/stiffCoeff(stiffOrder)));
set(handles.stiffOrderstr,'val',stiffOrder);
calibration=handles.calibrationORG;
xoffset=handles.xoffsetORG;
baseline=handles.baselineORG;
bL=handles.bL;
%stiffness=handles.stiffnessORG;
guidata(hObject, handles);
bp=handles.basep;
set(handles.BL01,'string',num2str(bL(1)));
set(handles.BL02,'string',num2str(bL(2)));
set(handles.bLwlc,'string',num2str(baseline));
set(handles.xOffset,'string',num2str(xoffset));
set(handles.cAlibration,'string',num2str(prod(calibration)));
set(handles.sTiffness,'string',num2str(1/stiffCoeff(stiffOrder)));
guidata(hObject, handles);

if(handles.FilterIndexData==1)
    
    Extension=handles.Extension;
    Return=handles.Return;
    
    
    Data=handles.RawData;
    
    tData(1,:)=Data(:,1);
    xData(1,:)=Data(:,2);
    FData(1,:)=Data(:,3);
    FData(2,:)=Data(:,6);
    
    
    [FAData,xAData]=AnalayzeData(bp,FData,xData,calibration,baseline,bL,stiffCoeff,xoffset);
    
    handles.FAData=FAData;
    handles.xAData=xAData;
    handles.tData=tData;
    
    xeData=xAData(Extension);
    xrData=xAData(Return);
    FeData(1,:)=FAData(Extension);
    FrData(1,:)=FAData(Return);
    
    handles.xrData=xrData;
    handles.FrData=FrData;
    handles.xeData=xeData;
    handles.FeData=FeData;
    handles.calibration=calibration;
    handles.xoffset=xoffset;
    handles.baseline=baseline;
    guidata(hObject, handles);
    axes(handles.DNA);
    cla(handles.DNA)
    plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
    plot(xrData,FrData,'or','MarkerSize',4)
    ymax=max(handles.FeData)+5;
    ymin=0;
    xmax=max(xeData)+0.02;
    xmax=max(xeData)+0.02;
    xmin=min(xeData)-0.01;
    axis([xmin xmax ymin ymax]);
    
    [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,handles.ZeroRows);
    handles.StrRetnData=StrRetnData;
    handles.tData_FFE=tData_FFE;
    handles.tData_FFR=tData_FFR;
    guidata(hObject, handles);
    AutoSaveVal = get(handles.AutoSave, 'Value');
    if(AutoSaveVal==1)
        [pathstr,DataName,ext] = fileparts(handles.filepartsData);
        StrRertnFileName=sprintf('%s-anlyzd.txt',DataName);
        cd(handles.PathNameData);
        dlmwrite(StrRertnFileName,StrRetnData,'delimiter','\t','newline','pc');
    end
    
    
    switch(handles.DataType)
        case 2
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFE(:,1),tData_FFE(:,2),tData_FFE(:,1),tData_FFE(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFE(:,3))+5],'YTick',[0:10:max(tData_FFE(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFE(:,2))-0.01 max(tData_FFE(:,2))+0.01],'YTick',[0:0.01:max(tData_FFE(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
            end
        case 3
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            if(AutoSaveVal==1)
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
            end
        case 4
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
            end
    end
    cd(handles.orgdir)
else
    
    xeDNA=handles.RawXeDNA;
    xrDNA=handles.RawXrDNA;
    feDNA=handles.RawFeDNA;
    frDNA=handles.RawFrDNA;
    
    
    
    [FrDNA,xrDNA]=AnalayzeData(bp,frDNA,xrDNA,calibration,baseline,bL,stiffCoeff,xoffset);
    [FeDNA,xeDNA]=AnalayzeData(bp,feDNA,xeDNA,calibration,baseline,bL,stiffCoeff,xoffset);
    
    
    cla(handles.DNA)
    axes(handles.DNA)
    plot(xrDNA,FrDNA,'ob','MarkerSize',4);
    plot(xeDNA,FeDNA,'ob','MarkerSize',4,'MarkerFaceColor','b');
    
    ymax=max(max(handles.FeDNA))+5;
    ymin=0;
    xmax=max(xeDNA)+0.02;
    xmin=min(xeDNA)-0.01;
    axis([xmin xmax ymin ymax]);
    set(handles.DNAtick,'value',1);
    [xrDNA i]=sort(xrDNA);
    FrDNA=FrDNA(i);
    
    handles.xeDNA=xeDNA;
    handles.xrDNA=xrDNA;
    handles.FeDNA=FeDNA;
    handles.FrDNA=FrDNA;
    guidata(hObject, handles);
    AutoSaveVal = get(handles.AutoSave, 'Value')
    
    
    if(AutoSaveVal==1)
        DNAonly=handles.FilePartsDNA;
        [pathstr,DNAname,ext] = fileparts(handles.FilePartsDNA);
        analyzedDNA= sprintf('%s-analyzed.txt',DNAname);
        M=reshapeData([xeDNA;FeDNA],[xrDNA;FrDNA]);
        cd(handles.PathName);
        dlmwrite(analyzedDNA,M,'delimiter','\t','newline','pc')
        set(handles.display,'string','Analyzed DNA Only File was saved');
        cd(handles.orgdir)
    end
    cd(handles.orgdir)
end




% --- Executes during object creation, after setting all properties.
function Legend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Legend.
function Legend_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function midpoint_Callback(hObject, eventdata, handles)

Salt_str=get(handles.midpoint, 'string');
SaltConc=str2double(Salt_str);
if(SaltConc>250 || SaltConc<2.6)
    warndlg('Salt concentration is out of acceptable range')
else
    Fmp=ForceMidpoint(SaltConc);
    FmpString=sprintf('Force Midpoint at %0.2f pN',Fmp);
    set(handles.display,'string',FmpString);
end
% hObject    handle to midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of midpoint as text
%        str2double(get(hObject,'String')) returns contents of midpoint as a double


% --- Executes on button press in Scaler.
function Scaler_Callback(hObject, eventdata, handles)
if(handles.DNAanalysedFlag==1 )
    xeDNA=handles.xeDNA;
    FeDNA=handles.FeDNA;
    minplat=handles.minplat;
    maxplat=handles.maxplat;
    Fmp=handles.Fmp;
    axes(handles.DNA)
    line([xeDNA(minplat),xeDNA(minplat)],[min(FeDNA),max(FeDNA)],'color','g','LineStyle',':','linewidth',3)
    line([xeDNA(maxplat),xeDNA(maxplat)],[min(FeDNA),max(FeDNA)],'color','g','LineStyle',':','linewidth',3)
    line([min(xeDNA),max(xeDNA)],[Fmp,Fmp],'color','m','LineStyle',':','linewidth',3)
    
    ymax=max(handles.FeDNA)+5;
    ymin=0;
    xmax=max(xeDNA)+0.02;
    xmin=min(xeDNA)-0.01;axis([xmin xmax ymin ymax]);
end

% --- Executes on button press in WLCbut.
function WLCbut_Callback(hObject, eventdata, handles)
if(handles.DNAanalysedFlag==1 )
    ds(1)= 0.3419;
    ds(2)= 0.0900;
    ds(3)=  1.3615e+03;
    
    WLC  = @(f) ds(1)*(1- 0.5*(ds(2)./f).^0.5 + (1/ds(3)).*f);
    fwlc=0.01:0.1:90;
    ewlc=WLC(fwlc);
    handles.wlcplot=plot(ewlc,fwlc,'r');
end
% hObject    handle to WLCbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WLCbut


% --- Executes on selection change in FFmenu.
function FFmenu_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if(handles.FilterIndexData==1)
    DataType=handles.DataType;
    switch val
        case 1
            if(DataType==2 || DataType==4)
                tData_FFE=handles.tData_FFE ;
                axes(handles.FFgraph)
                % delete(handles.H1);
                % delete(handles.H2);
                [AX,H1,H2]= plotyy(tData_FFE(:,1),tData_FFE(:,2),tData_FFE(:,1),tData_FFE(:,3));
                set(AX(2),'Ylim',[0 max(tData_FFE(:,3))+5],'YTick',[0:10:max(tData_FFE(:,3))+5],'YColor',[0.2,1,0.4]);
                set(AX(1),'box','off')
                set(AX(1),'Ylim',[min(tData_FFE(:,2))-0.01 max(tData_FFE(:,2))+0.01],'YTick',[0:0.01:max(tData_FFE(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
                ylabel(AX(2),'Force (pN)')
                ylabel(AX(1),'Extension (nm/bp)')
                xlabel('time (s)')
                handles.H1=H1;
                handles.H2=H2;
                guidata(hObject, handles);
            else
                axes(handles.FFgraph)
                delete(handles.H2);
                delete(handles.H1);
            end
            
        case 2
            
            if(DataType==3 || DataType==4)
                tData_FFR=handles.tData_FFR;
                axes(handles.FFgraph)
                % delete(handles.H2);
                %delete(handles.H1);
                
                [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
                set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5],'YColor',[0.2,1,0.4]);
                set(AX(1),'box','off')
                set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01],'YColor',[0,0.8,1],'XColor','w');
                ylabel(AX(2),'Force (pN)')
                ylabel(AX(1),'Extension (nm/bp)')
                xlabel('time (s)')
                handles.H1=H1;
                handles.H2=H2;
                guidata(hObject, handles);
                
                %plot
            else
                axes(handles.FFgraph)
                delete(handles.H2);
                delete(handles.H1);
            end
            
    end
end



% hObject    handle to FFmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FFmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FFmenu


% --- Executes during object creation, after setting all properties.
function FFmenu_CreateFcn(hObject, eventdata, handles)


% hObject    handle to FFmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------


% --- Executes on button press in SaveAs.
function SaveAs_Callback(hObject, eventdata, handles)
cd(handles.PathName)
if(handles.DataFileflag==0);
    [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save analyzed DNA Only file as');
    cd(handles.orgdir)
    if(FilterIndex==1)
        cd(handles.orgdir)
        xeDNA=handles.xeDNA;
        xrDNA=handles.xrDNA;
        FeDNA=handles.FeDNA;
        FrDNA=handles.FrDNA;
        M=reshapeData([xeDNA;FeDNA],[xrDNA;FrDNA]);
        cd(SaveFilePath);
        dlmwrite(FileName,M,'delimiter','\t','newline','pc')
        set(handles.display,'string','Analyzed DNA Only File was saved');
        cd(handles.orgdir)
        
    end
else
    tData=handles.tData;
    xAData=handles.xAData;
    FAData=handles.FAData;
    cd(handles.orgdir)
    [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,handles.ZeroRows);
    cd(handles.PathName)
    [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save analysed data as');
    
    if(FilterIndex==1)
        cd(SaveFilePath);
        dlmwrite(FileName,StrRetnData,'delimiter','\t','newline','pc');
        
        switch(handles.DataType)
            case 2
                cd(handles.PathName)
                [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save Force Feed back on Extension data as');
                if(FilterIndex==1)
                    cd(SaveFilePath);
                    dlmwrite(FileName,tData_FFE,'delimiter','\t','newline','pc');
                    cd(handles.orgdir)
                else
                    %timda data wasnt saved
                    cd(handles.orgdir)
                end
                
            case 3
                cd(handles.PathName)
                [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save Force Feed back on Return data as');
                if(FilterIndex==1)
                    cd(SaveFilePath);
                    dlmwrite(FileName,tData_FFR,'delimiter','\t','newline','pc');
                    cd(handles.orgdir)
                else
                    %timda data wasnt saved
                    cd(handles.orgdir)
                end
                
                
            case 4
                cd(handles.PathName)
                [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save Force Feed back on Extension data as');
                if(FilterIndex==1)
                    cd(SaveFilePath);
                    dlmwrite(FileName,tData_FFE,'delimiter','\t','newline','pc');
                else
                end
                
                [FileName,SaveFilePath,FilterIndex] = uiputfile('.txt','Save Force Feed back on Return data as');
                if(FilterIndex==1)
                    cd(SaveFilePath);
                    dlmwrite(FileName,tData_FFR,'delimiter','\t','newline','pc');
                else
                    %timda data wasnt saved
                    cd(handles.orgdir)
                end
        end
    else
        %Data wasnt saved
        cd(handles.orgdir)
    end
    cd(handles.orgdir)
end



function pullRate_Callback(hObject, eventdata, handles)


% hObject    handle to pullRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pullRate as text
%        str2double(get(hObject,'String')) returns contents of pullRate as a double


% --- Executes during object creation, after setting all properties.
function pullRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pullRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stiffOrderstr.
function stiffOrderstr_Callback(hObject, eventdata, handles)
if(handles.DNAanalysedFlag==1)
    calibration_str=get(handles.cAlibration,'string');
    calibration=str2double(calibration_str);
    stiffFile=handles.stiffFile;
    axes(handles.DNA);
    cla(handles.DNA)
    n=find(all(stiffFile==0,2));
    Stext(1,:)=stiffFile((1:n-1),2);
    StFext(1,:)=stiffFile((1:n-1),3);
    StFext(2,:)=stiffFile((1:n-1),6);
    StFext=calibration.*mean(StFext);
    stiffOrder=get(handles.stiffOrderstr,'val');
    p= findStiffnessFit(calibration,stiffFile,stiffOrder);
    XX=correctStiffness(p,Stext,StFext,1000,1);
    set(handles.sTiffness,'string',num2str(1/p(stiffOrder)))
    plot(XX,StFext,'ob','MarkerSize',4,'MarkerFaceColor','b');hold on
    ymax=max(StFext)+5;
    ymin=min(StFext)-5;
    xmax=max(XX)+1;
    xmin=min(XX)-1;
    axis([xmin xmax ymin ymax]);
    handles.stiffCoeff=p;
    guidata(hObject, handles);
end
% hObject    handle to stiffOrderstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stiffOrderstr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stiffOrderstr


% --- Executes during object creation, after setting all properties.
function stiffOrderstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stiffOrderstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
winopen('help.pdf')
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipushtool12_ClickedCallback(hObject, eventdata, handles)
winopen('help.pdf')
% hObject    handle to uipushtool12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sTiffness_Callback(hObject, eventdata, handles)

% hObject    handle to sTiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sTiffness as text
%        str2double(get(hObject,'String')) returns contents of sTiffness as a double


% --- Executes on button press in stiffButton.
function stiffButton_Callback(hObject, eventdata, handles)
if(handles.stiffflag==1)
    if(handles.DNAanalysedFlag==1 )
        p=handles.stiffCoeff;
        stifftext = StiffnessText(p);
        set(handles.display,'string',stifftext,'fontsize',8);
        
    else
        warndlg('Analyze the DNA first')
    end
else
    warndlg('Load the stiffness file')
end
% hObject    handle to stiffButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Overlap.
function Overlap_Callback(hObject, eventdata, handles)
%stiffOrder=handles.stiffOrder;

if(handles.FilterIndexData==1)
    stiffCoeff=handles.stiffCoeff;
    calibration=handles.calibration;
    xoffset=handles.xoffset
    baseline=handles.baseline;
    bL=handles.bL;
    bp=handles.basep;
    Data=handles.RawData;
    Extension=handles.Extension;
    Return=handles.Return;
    tData(1,:)=Data(:,1);
    xData(1,:)=Data(:,2);
    FData(1,:)=Data(:,3);
    FData(2,:)=Data(:,6);
    overlapVal=get(handles.overlapmenu,'val');
    
    for(i=1:5)
        [FAData,xAData]=AnalayzeData(bp,FData,xData,calibration,baseline,bL,stiffCoeff,xoffset);
        [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,handles.ZeroRows);
        [df,dx] = OverLap(StrRetnData(:,2),StrRetnData(:,1),overlapVal);
        xoffset=xoffset+dx;
        baseline=baseline+df;
        xAData = xAData+dx;
        FAData = FAData+df;
    end
    handles.StrRetnData=StrRetnData;
    handles.tData_FFE=tData_FFE;
    handles.tData_FFR=tData_FFR;
    set(handles.bLwlc,'string',num2str(baseline));
    set(handles.xOffset,'string',num2str(xoffset));
    handles.tData = tData;
    handles.xAData= xAData;
    handles.FAData= FAData;
    
    xeData=xAData(Extension);
    xrData=xAData(Return);
    FeData(1,:)=FAData(Extension);
    FrData(1,:)=FAData(Return);
    
    handles.xrData=xrData;
    handles.FrData=FrData;
    handles.xeData=xeData;
    handles.FeData=FeData;
    handles.xoffset=xoffset;
    handles.baseline=baseline;
    
    guidata(hObject, handles);
    axes(handles.DNA);
    cla(handles.DNA) ;
    plot(xeData,FeData,'or','MarkerSize',4,'MarkerFaceColor','r');hold on
    plot(xrData,FrData,'or','MarkerSize',4)
    ymax=max(handles.FeData)+5;
    ymin=0;
    
    xmax=max(xeData)+0.02;
    xmin=min(xeData)-0.01;
    axis([xmin xmax ymin ymax]);
    [StrRetnData,tData_FFE,tData_FFR]=SaveData(tData,xAData,FAData,handles.ZeroRows);
    AutoSaveVal = get(handles.AutoSave, 'Value');
    if(AutoSaveVal==1)
        [pathstr,DataName,ext] = fileparts(handles.filepartsData);
        StrRertnFileName=sprintf('%s-anlyzd.txt',DataName);
        cd(handles.PathNameData);
        dlmwrite(StrRertnFileName,StrRetnData,'delimiter','\t','newline','pc');
    end
    
    
    switch(handles.DataType)
        case 2
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFE(:,1),tData_FFE(:,2),tData_FFE(:,1),tData_FFE(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFE(:,3))+5],'YTick',[0:10:max(tData_FFE(:,3))+5],'YColor',[0.2,1,0.4]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFE(:,2))-0.01 max(tData_FFE(:,2))+0.01],'YTick',[0:0.01:max(tData_FFE(:,2))+0.01]);
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                set(handles.display,'string','Analyzed Data File was updated');
            end
        case 3
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01]);
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            if(AutoSaveVal==1)
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
                set(handles.display,'string','Analyzed Data File was updated');
            end
        case 4
            delete(handles.H2);
            delete(handles.H1);
            axes(handles.FFgraph)
            [AX,H1,H2]= plotyy(tData_FFR(:,1),tData_FFR(:,2),tData_FFR(:,1),tData_FFR(:,3));
            set(AX(2),'Ylim',[0 max(tData_FFR(:,3))+5],'YTick',[0:10:max(tData_FFR(:,3))+5]);
            set(AX(1),'box','off')
            set(AX(1),'Ylim',[min(tData_FFR(:,2))-0.01 max(tData_FFR(:,2))+0.01],'YTick',[0:0.01:max(tData_FFR(:,2))+0.01]);
            ylabel(AX(2),'Force (pN)')
            ylabel(AX(1),'Extension (nm/bp)')
            xlabel('time (s)')
            handles.H1=H1;
            handles.H2=H2;
            guidata(hObject, handles);
            
            if(AutoSaveVal==1)
                FfeFileName=sprintf('%s-anlyzd-FFEtime.txt',DataName);
                handles.FfeFileName=FfeFileName;
                dlmwrite(FfeFileName,tData_FFE,'delimiter','\t','newline','pc');
                FfrFileName=sprintf('%s-anlyzd-FFRtime.txt',DataName);
                dlmwrite(FfrFileName,tData_FFR,'delimiter','\t','newline','pc');
                set(handles.display,'string','Analyzed Data File was updated');
            end
    end
    cd(handles.orgdir)
    
end

guidata(hObject, handles);

% hObject    handle to Overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function bLwlc_Callback(hObject, eventdata, handles)
% hObject    handle to bLwlc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bLwlc as text
%        str2double(get(hObject,'String')) returns contents of bLwlc as a double


% --- Executes during object creation, after setting all properties.
function bLwlc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bLwlc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function midpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to midpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function BL01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BL01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function BL02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BL02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BL01_Callback(hObject, eventdata, handles)
% hObject    handle to BL01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BL01 as text
%        str2double(get(hObject,'String')) returns contents of BL01 as a double



function BL02_Callback(hObject, eventdata, handles)
% hObject    handle to BL02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BL02 as text
%        str2double(get(hObject,'String')) returns contents of BL02 as a double



function xOffset_Callback(hObject, eventdata, handles)
% hObject    handle to xOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xOffset as text
%        str2double(get(hObject,'String')) returns contents of xOffset as a double


% --- Executes during object creation, after setting all properties.
function xOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cAlibration_Callback(hObject, eventdata, handles)
% hObject    handle to cAlibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cAlibration as text
%        str2double(get(hObject,'String')) returns contents of cAlibration as a double


% --- Executes during object creation, after setting all properties.
function cAlibration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cAlibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function sTiffness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sTiffness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function stiffButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stiffButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.


% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function AnaDNA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnaDNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function bp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function tolerance_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as a double



function bp_Callback(hObject, eventdata, handles)
% hObject    handle to bp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bp as text
%        str2double(get(hObject,'String')) returns contents of bp as a double


% --- Executes on button press in AutoSave.
function AutoSave_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutoSave


% --- Executes during object creation, after setting all properties.
function AutoSave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in overlapmenu.
function overlapmenu_Callback(hObject, eventdata, handles)
% hObject    handle to overlapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns overlapmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from overlapmenu


% --- Executes during object creation, after setting all properties.
function overlapmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pN_Callback(hObject, eventdata, handles)
% hObject    handle to pN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pN as text
%        str2double(get(hObject,'String')) returns contents of pN as a double
if(handles.DNAanalysedFlag==1 )
    Micrn=handles.feDNARaw;
    pNs=handles.FeDNA;
    pNVal_str=get(handles.pN,'string');
    pNVal=str2double(pNVal_str);
    
    if( pNVal>max(pNs)|| pNVal < min(pNs))
        warndlg('Value entered is out of bound')
    else
        
        pNInd=find(pNs>=pNVal);
        y2=Micrn(pNInd(1));
        y1=Micrn(pNInd(1)-1);
        x2=pNs(pNInd(1));
        x1=pNs(pNInd(1)-1);
        
        MicrnVal=pNVal*(y1+y2)/(x1+x2);
        set(handles.Mic,'string',num2str(MicrnVal));
    end
else
    warndlg('Please analyse the DNA first')
    
end

% --- Executes during object creation, after setting all properties.
function pN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mic_Callback(hObject, eventdata, handles)
% hObject    handle to Mic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mic as text
%        str2double(get(hObject,'String')) returns contents of Mic as a double


% --- Executes during object creation, after setting all properties.
function Mic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

ha = axes('units','normalized', 'position',[0 0 1 1]);
% Move the background axes to the bottom

uistack(ha,'bottom');
I=imread('0DNVzD2.jpg');
imagesc(I)

% --- Executes on button press in SplitAndSave.
function SplitAndSave_Callback(hObject, eventdata, handles)
% hObject    handle to SplitAndSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd(handles.PathNameData)
[FF,stpath,ffFlag] = uigetfile({'*.*','All Files(*.*)'}, 'Chose the Force Feedback file ');
[pathstr,DataName,ext] = fileparts(handles.filepartsData);
b = strcat(stpath,FF);
timedata=importdata(b);
t=timedata(:,1);
x=timedata(:,2);
f=timedata(:,3);
cd(handles.orgdir)
SeperateForce(t,x,f,100,handles.PathNameData,DataName,handles.orgdir)

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

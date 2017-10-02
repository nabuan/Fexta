function [ImExtPx] = TrackExtPx(Bead_1,Bead_2,Bead_1R,Bead_2R,IntThreshold,datpath,AxesHandle,MotionRange)

Tracker=AxesHandle;

datpath=sprintf('%s\\',datpath);
filenm=dir([datpath,'*.jpg']);
filenm={filenm(:).name};



Cthr=IntThreshold;
Cthr=0.8;


I_t=[];
I_R=MotionRange;
xR=I_R(1);dxR=I_R(2);yR=I_R(3);dyR=I_R(4);
xL=Bead_1R(1);dxL=Bead_1R(2);yL=Bead_1R(3);dyL=Bead_1R(4); 
xS=Bead_2R(1);dxS=Bead_2R(2);yS=Bead_2R(3),dyS=Bead_2R(4); % bead2 cordinates

h_im=findobj(findobj(Tracker,'Type','Image'));
hold(Tracker,'on');
if isempty(h_im)
    imshow(I_t,[],'Parent',Tracker);
    h_im=findobj(findobj(Tracker,'Type','Image'));
else
    h_im.CData=I_t;
end

ImExtPx=nan(numel(filenm),4);
offset_small=0.5.*(size(Bead_2));
offset_large=0.5.*(size(Bead_1));


%stack all images and correct bacground and let the corrected be I_t_crcted
% find Bead_2 ref image from the corrected images
%Bead_2=I_t_crcted(xR:xR+dxR,yR:yR+dyR)




for i=1:numel(filenm)
    I_t=double(imread([datpath,filenm{i}]));
    I_t=I_t(xR:xR+dxR,yR:yR+dyR);
    
%     if(i>1)
%     I_prev=double(imread([datpath,filenm{i-1}]));
%     Bead_1=I_prev(Bead_1px-offset_large(2):Bead_1px-offset_large(2)+dxL,row_large-offset_large(1):row_large-offset_large(1)+dyL);
%     Bead_2=I_prev(xS:xS+dxS,yS:yS+dyS);
%     end
    
    
    C_large=normxcorr2(Bead_1,I_t);
    CthrL=0.8*max(max(C_large));
    [is,js,vs]=find((C_large>CthrL));
    row_large=sum(is.*vs)./sum(vs);
    Bead_1px=sum(js.*vs)./sum(vs);
    
    C_small=normxcorr2(Bead_2,I_t);
    CthrS=0.8*max(max(C_small));
    [is,js,vs]=find((C_small> CthrS));
    row_small=sum(is.*vs)./sum(vs);
    Bead_2px=sum(js.*vs)./sum(vs);
    
    [~,name,ext] = fileparts(filenm{i});
    FileInd =str2double(strtok(name, 'Img'));
    
    ImExtPx(i,1)=FileInd;
    ImExtPx(i,2)= Bead_1px;
    ImExtPx(i,3)= Bead_2px;
    ImExtPx(i,4)= abs(Bead_2px - Bead_1px) ;
    h_im.CData=I_t;
    drawnow;
    
    h_large=findobj(Tracker.Children,'Tag','RedCross');
    h_small=findobj(Tracker.Children,'Tag','BlueCross');
    if isempty(h_large) || isempty(h_small)
        delete(findobj(Tracker,'Type','line'))
        h_large=plot(Bead_1px-offset_large(2),row_large-offset_large(1),'xr','MarkerSize',15,'Parent',Tracker,'Tag','RedCross');hold on;drawnow;
        h_small=plot(Bead_2px-offset_small(2),row_small-offset_small(1),'xb','MarkerSize',15,'Parent',Tracker,'Tag','BlueCross');hold on;drawnow;
    else
        h_large.XData=Bead_1px-offset_large(2);h_large.YData=row_large-offset_large(1);
        h_small.XData=Bead_2px-offset_small(2);h_small.YData=row_small-offset_small(1);
    end
    
end
end


function [ImExtPx] = TrackExtPx(Bead_1,Bead_2,IntThreshold,datpath,AxesHandle,MotionRange)

Tracker=AxesHandle;
%datpath=uigetdir;
datpath=sprintf('%s\\',datpath);
filenm=dir([datpath,'*.jpg']);
filenm={filenm(:).name};
%I_ref=double(imread([datpath,filenm{1}]));

%figure(1)
%imshow(I_ref,[]);
%title('Chose the range of motion')
Cthr=IntThreshold;
Cthr=0.8;

% rectR= getrect(figure(1));
% rectR=ceil(rectR);yR=rectR(1);xR=rectR(2);dyR=rectR(3);dxR=rectR(4);
% I_R=I_ref(xR:xR+dxR,yR:yR+dyR);
% close(figure(1));

I_R=MotionRange;
xR=I_R(1);dxR=I_R(2);yR=I_R(3);dyR=I_R(4);

h_im=findobj(findobj(Tracker,'Type','Image'));
hold(Tracker,'on');
if isempty(h_im)
    imshow(I_R,[],'Parent',Tracker);
    h_im=findobj(findobj(Tracker,'Type','Image'));
else
    h_im.CData=I_R;
end

ImExtPx=nan(numel(filenm),4);
offset_small=0.5.*(size(Bead_2));
offset_large=0.5.*(size(Bead_1));
% hf=figure(1);
% h_1=subplot(3,2,1);hold on
% h_2=subplot(3,2,4);
% h_2.XLim=[0 numel(filenm)];
% h_2.YLim=[0 500*scale];axis square;box on;drawnow;hold on
% h_2.XLabel.String='Step';
% h_2.YLabel.String='Distance';
% h_3=subplot(3,2,2);
% h_3.YLim=[0 30*scale];
% h_3.XLim=[100*scale 220*scale];axis square;box on;drawnow;hold on
% h_3.XLabel.String='Extension';
% h_3.YLabel.String='Force';
% h_4=subplot(3,2,3);hold on
% h_5=subplot(3,2,5);hold on
for i=1:numel(filenm)
    I_t=double(imread([datpath,filenm{i}]));
    I_t=I_t(xR:xR+dxR,yR:yR+dyR);
    
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
    %     if isempty(h_im)
    %         h_im=imshow(I_t,[],'Parent',h_1);hold on
    %         h_im.Tag='image';
    %         h_large=plot(col_large(i)-offset_large(2),row_large-offset_large(1),'xr','MarkerSize',15,'Parent',h_1);hold on;drawnow;
    %         h_small=plot(col_small(i)-offset_small(2),row_small-offset_small(1),'xb','MarkerSize',15,'Parent',h_1);hold on;drawnow;
    %         h_imCL=imshow(C_large,[],'Parent',h_4);hold on
    %         h_imSM=imshow(C_small,[],'Parent',h_5);hold on
    %     else
    %         h_im.CData=I_t;
    %         h_imCL.CData=C_large;
    %         h_imSM.CData=C_small;
    %         h_large.XData=col_large(i)-offset_large(2);h_large.YData=row_large-offset_large(1);
    %         h_small.XData=col_small(i)-offset_small(2);h_small.YData=row_small-offset_small(1);
    %     end
    %     drawnow;
    %
    %     dist_large(i)=abs(col_large(i)-col_large(1));
    %     dist(i)=abs(col_large(i)-col_small(i));
    %     plot(i,dist_large(i),'or','Parent',h_2),hold on;
    %     plot(i, dist(i),'oc','Parent',h_2),hold on;
    %     m=redbluecmap(2);
    %     m(2,:)=[];
    %     t=str2num(filenm{i}(4));
    %     mk={'x','o'}
    %     plot(dist(i),dist_large(i),'Parent',h_3,'Color',m(t,:),'Marker',mk{t}),hold on;
    %     drawnow;
    %     %pause(0.5);
end
end


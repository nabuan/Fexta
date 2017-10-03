function [ImExtPx] = TrackExtPx(Bead_1,Bead_2,Bead_1R,Bead_2R,IntThreshold,datpath,AxesHandle,MotionRange)


filenm=dir([datpath,'*.jpg']);
filenm={filenm(:).name};

%Renaming variables for internal use:
Tracker=AxesHandle;
I_R=MotionRange;
xR=I_R(1);dxR=I_R(2);yR=I_R(3);dyR=I_R(4);
xL=Bead_1R(1);dxL=Bead_1R(2);yL=Bead_1R(3);dyL=Bead_1R(4);
xS=Bead_2R(1);dxS=Bead_2R(2);yS=Bead_2R(3);dyS=Bead_2R(4);

ImExtPx=nan(numel(filenm),4);
offset_small=0.5.*(size(Bead_2));
offset_large=0.5.*(size(Bead_1));

%{
if isplot
    h_im=findobj(findobj(Tracker,'Type','Image'));
    hold(Tracker,'on');
    if isempty(h_im)
        imshow(I_t,[],'Parent',Tracker);
        h_im=findobj(findobj(Tracker,'Type','Image'));
    else
        h_im.CData=I_t;
    end
end
%}

%Loading sequence into a stack
dirname='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/SameTypeBead/01DNA/dsDNA/';
%dirname='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/Images/';
fname=dir([dirname,'*.jpg']);
IM=imread([dirname,fname(1).name]);
IMstack=zeros([size(IM),numel(fname)]);
IMstack(:,:,1)=IM;
for f=2:numel(fname)
    IM=imread([dirname,fname(f).name]);
    IMstack(:,:,f)=IM;
end
%Reduce size by 4 fold
%Remove mean across time from each pixel
%and perform simple gaussian filtering over the result

redfact=4;
IMred=Reduce_Stack(IMstack,redfact,redfact,1,'Mean');

for f=1:size(IMred,3)
    IMred(:,:,f) = imgaussfilt(IMred(:,:,f),2);
end

IMred2=IMred;
IMred=bsxfun(@minus,IMred,mean(IMred,3));

BGcorr=bsxfun(@minus,IMred,mean(IMred,3));
%Click on approximate center of the bead
figure(999);clf(999);
imshow(BGcorr(:,:,1),[]);
im_pad=15;search_pad=25;%search_pad should be larger than im_pad
tb_midpt_xy=round(ginput(1));
tb_im=BGcorr((tb_midpt_xy(2)-search_pad):(tb_midpt_xy(2)+search_pad),(tb_midpt_xy(1)-search_pad):(tb_midpt_xy(1)+search_pad),1);
close(999)

%Readjust reference image
[tb_midpt_relxy,~,~]=imfindcircles(tb_im,[7 15],'ObjectPolarity','dark');
tb_midpt_xy=round((tb_midpt_relxy-1)+(tb_midpt_xy-search_pad));
tb_im=IMred((tb_midpt_xy(2)-im_pad):(tb_midpt_xy(2)+im_pad),(tb_midpt_xy(1)-im_pad):(tb_midpt_xy(1)+im_pad),1);

figure(1000);clf(1000);
imagesc(tb_im);colorbar;
colormap('gray');axis equal;caxis([min(tb_im(:)) max(tb_im(:))])
drawnow;

loc=nan(size(IMred,3),2);
ref=nan(2*im_pad+1,2*im_pad+1,size(IMred,3));
ref(:,:,1)=tb_im;
for f=2:size(IMred,3)
    srch=IMred(...
        (tb_midpt_xy(2)-search_pad):(tb_midpt_xy(2)+search_pad),...
        (tb_midpt_xy(1)-search_pad):(tb_midpt_xy(1)+search_pad),...
        f);
    cc=normxcorr2(ref(:,:,f-1),srch);
    
    figure(1);clf(1)
    imagesc(srch);colorbar;colormap('gray');axis equal;caxis([min(srch(:)) max(srch(:))])
    drawnow;
    
    %Portions without complete overlap of ref and srch are set to 0
    msk=true(size(cc));msk(size(ref,1):size(ref,1)+(size(srch)-size(ref,1)),...
        size(ref,2):size(ref,2)+(size(srch)-size(ref,2)))=false;
    cc(msk)=0;cc(cc<=0.8*nanmax(cc(:)))=0;
    
    [is,js,vs]=find(cc);
    cc_midpt_relxy=[sum(js.*vs)./sum(vs),sum(is.*vs)./sum(vs)];
    
    %midpoint of matrix cc is centered at tb_midpt_xy
    cc_midpt_xy=[cc_midpt_relxy(1)-round(size(cc,1)./2),...
        cc_midpt_relxy(2)-round(size(cc,2)./2)];
    
    loc(f,:)=cc_midpt_xy+tb_midpt_xy;
    tb_midpt_xy=round(loc(f,:));
    if mod(f,10)==0
        ref(:,:,f)=IMred(...
            (tb_midpt_xy(2)-im_pad):(tb_midpt_xy(2)+im_pad),...
            (tb_midpt_xy(1)-im_pad):(tb_midpt_xy(1)+im_pad),...
            f);
    else
        ref(:,:,f)=ref(:,:,f-1);
    end
    figure(2);clf(2)
    temp=ref(:,:,f);
    imagesc(ref(:,:,f));colorbar;colormap('gray');axis equal;caxis([min(temp(:)) max(temp(:))])
    drawnow;
end

% im3dscroll(ref,[],[1:size(ref,3)]);colorbar;
% colormap('gray');axis equal;caxis([-50 50])
% 
% %Check on the image:
% im3dscroll(IMred,[loc(:,[2,1]),[1:size(IMred,3)]'],[1:size(IMred,3)]);colorbar;
% colormap('gray');axis equal;caxis([-50 50])

%Check on the image:
im3dscroll(IMstack,[4*loc(:,[2,1]),[1:size(IMred,3)]'],[1:size(IMred,3)]);colorbar;
colormap('gray');axis equal;caxis([0 300])


%}
%{

im3dscroll(IMstack,[],1:size(IMstack,3));

for i=1:numel(filenm)
    I_t=double(imread([datpath,filenm{i}]));
    I_t=I_t(xR:xR+dxR,yR:yR+dyR);
    
%{
    if(i>1)
        I_prev=double(imread([datpath,filenm{i-1}]));
        Bead_1=I_prev(Bead_1px-offset_large(2):Bead_1px-offset_large(2)+dxL,row_large-offset_large(1):row_large-offset_large(1)+dyL);
        Bead_2=I_prev(xS:xS+dxS,yS:yS+dyS);
    end
%}
    
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
    
    ImExtPx(i,1)= FileInd;
    ImExtPx(i,2)= Bead_1px;
    ImExtPx(i,3)= Bead_2px;
    ImExtPx(i,4)= abs(Bead_2px - Bead_1px) ;
    
    if isplot
        drawnow;
        h_im.CData=I_t;
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
%end

%}
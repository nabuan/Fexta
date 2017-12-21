function [pos1,pos2]=TrackBead(dirname)
%This function takes directory name of the image frames as input and
%outputs the stage bead position in successive frames. Tracking is possible
%in images with with non-uniform illumination. The issue of bias in bead
%position due to accumulation of rounding and tracking errors is fixed in
%this version

%Algorithm summary:
%1. Calculate reference bead center using Hough transform
%2. Hough transform is performed on corrected image.
%3. Calculate corrected image by subtracting mean across frames from each pixel value in current frame
%4. Use normalized cross correlation using original/reduced image to find bead position in successive frames

%Algorithm 2
%1. Perform high-pass fourier filtering on 1D intensity profiles
%2. Detect peaks within expected

%Parameters:
updateref=1000;         %Steps after which reference image is recalculated
im_pad=30;              %Smallest value that can accomodate entire reference bead
search_pad=45;          %im_pad < search_pad. Smalles value to include bead position in next frame
redfact=2;              %Factor by which image is reduced. Above parameters should be rescaled automatically in future version.
isplot=false;
%--------------------------------------------------------------------------
%Load sequence into a stack
if isempty(dirname)
    dirname='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/SameTypeBead/01DNA/dsDNA/';
end
fname=dir([dirname,'*.jpg']);
IM=imread([dirname,fname(1).name]);
IMstack=zeros([size(IM),numel(fname)]);
fname={fname(:).name}';
k = strfind(fname{1},'.');
seq=zeros(numel(fname),1);
for f=1:numel(fname)
    seq(f)=str2double(fname{f}(k-4:k-1));
end

%Sequence of images should be 1xxxx, 3xxxx, 4xxxx and 2xxxx.
%xxxx reflects the correct sequence order.
[~,sorted]=sort(seq);
for f=1:numel(fname)
    IM=imread([dirname,fname{sorted(f)}]);
    IMstack(:,:,f)=IM;
end

IM_reduced=Reduce_Stack(IMstack,redfact,redfact,1,'Mean');
for f=1:size(IM_reduced,3)
    IM_reduced(:,:,f) = imgaussfilt(IM_reduced(:,:,f),2);
end
indmax=max([50,size(IM_reduced,3)]);
IM_corr=bsxfun(@minus,IM_reduced,mean(IM_reduced(:,:,1:indmax),3));

%--------------------------------------------------------------------------
%------------------------------STAGE BEAD----------------------------------
%--------------------------------------------------------------------------
%Click on approximate center of the stage bead
figure(999);clf(999);
imshow(IM_corr(:,:,1),[]);
bead_midpt_xy=round(ginput(1));
bead_im=IM_corr((bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),(bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),1);
close(999)

%Re-adjust reference image to find bead
[bead_midpt_relxy,~,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
[~,keep]=min(sum((bead_midpt_relxy-search_pad).^2,2));
bead_midpt_relxy=bead_midpt_relxy(keep,:);
bead_midpt_xy=round((bead_midpt_relxy-1)+(bead_midpt_xy-search_pad));
bead_im=IM_reduced((bead_midpt_xy(2)-im_pad):(bead_midpt_xy(2)+im_pad),(bead_midpt_xy(1)-im_pad):(bead_midpt_xy(1)+im_pad),1);

loc=nan(size(IM_reduced,3),2);
loc(1,:)=bead_midpt_xy;
ref=nan(2*im_pad+1,2*im_pad+1,size(IM_reduced,3));
ref(:,:,1)=bead_im;
round_err=[0 0];

for f=2:size(IM_reduced,3)
    srch=IM_reduced(...
        (bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),...
        (bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),...
        f);
    cc=normxcorr2(ref(:,:,f-1),srch);
    
    if isplot
        %Plot image in which bead is searched for
        figure(1);clf(1);movegui(1,'northeast')
        imagesc(srch);colorbar;colormap('gray');
        axis equal;caxis([min(srch(:)) max(srch(:))])
        drawnow;
        
        %Plot reference bead
        figure(2);clf(2);movegui(2,'southeast')
        refim=ref(:,:,f-1);axis equal;
        caxis([min(refim(:)) max(refim(:))])
        drawnow;
    end
    
    %Portions without complete overlap of ref and srch are set to 0
    msk=true(size(cc));msk(size(ref,1):size(ref,1)+(size(srch)-size(ref,1)),...
        size(ref,2):size(ref,2)+(size(srch)-size(ref,2)))=false;
    cc(msk)=0;cc(cc<=0.9*nanmax(cc(:)))=0;
    
    [is,js,vs]=find(cc);
    cc_midpt_relxy=[sum(js.*vs)./sum(vs),sum(is.*vs)./sum(vs)];
    
    %midpoint of matrix cc is centered at tb_midpt_xy
    cc_midpt_xy=[cc_midpt_relxy(1)-round(size(cc,1)./2),...
        cc_midpt_relxy(2)-round(size(cc,2)./2)];
    
    loc(f,:)=cc_midpt_xy+bead_midpt_xy+round_err;% Add rounding error from previous iteration
    round_err=loc(f,:)-round(loc(f,:));%This error is added in the next iteration
    bead_midpt_xy=round(loc(f,:));
    if mod(f,updateref)==0
        bead_im=IM_corr((bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),(bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),f);
        [bead_midpt_relxy,~,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
        bead_midpt_xy=round((bead_midpt_relxy-1)+(bead_midpt_xy-search_pad));
        
        ref(:,:,f)=IM_reduced(...
            (bead_midpt_xy(2)-im_pad):(bead_midpt_xy(2)+im_pad),...
            (bead_midpt_xy(1)-im_pad):(bead_midpt_xy(1)+im_pad),...
            f);
    else
        ref(:,:,f)=ref(:,:,f-1);
    end
end
loc=loc(:,[2,1]);
pos=redfact*loc;

%{
%Plot reference beads in all frames:
im3dscroll(ref,[],(0:size(ref,3)-1));colorbar;
colormap('gray');axis equal;caxis([0 300])

%Plot reduced image with tracked position: less memory/faster scrolling
im3dscroll(IM_reduced,[loc,(1:size(IM_reduced,3))'],(1:size(IM_reduced,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])

%Plot on original stack: more memory/slower scrolling
im3dscroll(IMstack,[pos,(1:size(IMstack,3))'],(0:size(IMstack,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])
%}

%pos contains Bead2 position in the first column
temp=pos;
temp(sorted,:)=pos;
pos=temp;

%Below vars only used for plotting plot
loc1=loc;
pos1=temp;
%--------------------------------------------------------------------------
%------------------------------TRAP BEAD-----------------------------------
%--------------------------------------------------------------------------
updatered=inf;

%Click on approximate center of the stage bead
im3dscroll(IM_reduced,[],(1:size(IM_reduced,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])
drawnow;hf=gcf;

keypressed=0;
while keypressed~=1
    keypressed=waitforbuttonpress;
end
h=gca;
IM_currplane=h.Children.CData;
bead_midpt_xy=round(ginput(1));
bead_im=IM_currplane((bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),(bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),1);
close(hf);

%Re-adjust reference image to find bead
[bead_midpt_relxy,r_estimated,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
bead_midpt_xy=round((bead_midpt_relxy-1)+(bead_midpt_xy-search_pad));
bead_im=IM_currplane((bead_midpt_xy(2)-im_pad):(bead_midpt_xy(2)+im_pad),(bead_midpt_xy(1)-im_pad):(bead_midpt_xy(1)+im_pad),1);

%Find estimate of bead boundary
bead_im_bin=false(size(bead_im));
bead_im_bin(im_pad+1,im_pad+1)=true;
bead_im_dist=bwdist(bead_im_bin);
bead_im_corrected=bead_im;

%Shuffling pixels to remove effect of beam/shadow
vals=bead_im_corrected(bead_im_dist<r_estimated*0.7);
bead_im_corrected(bead_im_dist<r_estimated*0.7)=vals(randperm(numel(vals)));
vals=bead_im_corrected(bead_im_dist>r_estimated*1.2);
bead_im_corrected(bead_im_dist>r_estimated*1.2)=vals(randperm(numel(vals)));

%Updating bead_im
bead_im=bead_im_corrected;



%--------------------End of fix--------------------------------------------
loc=nan(size(IM_reduced,3),2);
%Find loc in the first image of the set
f=1;round_err=[0 0];
srch=IM_reduced(...
    (bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),...
    (bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),...
    f);
cc=normxcorr2(bead_im,srch);

%Portions without complete overlap of ref and srch are set to 0
msk=true(size(cc));msk(size(ref,1):size(ref,1)+(size(srch)-size(ref,1)),...
    size(ref,2):size(ref,2)+(size(srch)-size(ref,2)))=false;
cc(msk)=0;cc(cc<=0.9*nanmax(cc(:)))=0;

[is,js,vs]=find(cc);
cc_midpt_relxy=[sum(js.*vs)./sum(vs),sum(is.*vs)./sum(vs)];
%midpoint of matrix cc is centered at bead_midpt_xy
cc_midpt_xy=[cc_midpt_relxy(1)-round(size(cc,1)./2),...
    cc_midpt_relxy(2)-round(size(cc,2)./2)];
loc(f,:)=cc_midpt_xy+bead_midpt_xy+round_err;% Add rounding error from previous iteration
bead_midpt_xy=round(loc(f,:));
%--------------------End of fix--------------------------------------------

ref=nan(2*im_pad+1,2*im_pad+1,size(IM_reduced,3));
ref(:,:,1)=bead_im;
round_err=[0 0];

for f=2:size(IM_reduced,3)
    srch=IM_reduced(...
        (bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),...
        (bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),...
        f);
    cc=normxcorr2(ref(:,:,f-1),srch);
    
    if isplot
        %Plot image in which bead is searched for
        figure(1);clf(1);movegui(1,'northeast')
        imagesc(srch);colorbar;colormap('gray');
        axis equal;caxis([min(srch(:)) max(srch(:))])
        drawnow;
        
        %Plot reference bead
        figure(2);clf(2);movegui(2,'southeast')
        refim=ref(:,:,f-1);axis equal;
        caxis([min(refim(:)) max(refim(:))])
        drawnow;
    end
    
    %Portions without complete overlap of ref and srch are set to 0
    msk=true(size(cc));msk(size(ref,1):size(ref,1)+(size(srch)-size(ref,1)),...
        size(ref,2):size(ref,2)+(size(srch)-size(ref,2)))=false;
    cc(msk)=0;cc(cc<=0.9*nanmax(cc(:)))=0;
    
    [is,js,vs]=find(cc);
    cc_midpt_relxy=[sum(js.*vs)./sum(vs),sum(is.*vs)./sum(vs)];
    
    %midpoint of matrix cc is centered at bead_midpt_xy
    cc_midpt_xy=[cc_midpt_relxy(1)-round(size(cc,1)./2),...
        cc_midpt_relxy(2)-round(size(cc,2)./2)];
    
    loc(f,:)=cc_midpt_xy+bead_midpt_xy+round_err;% Add rounding error from previous iteration
    round_err=loc(f,:)-round(loc(f,:));%This error is added in the next iteration
    bead_midpt_xy=round(loc(f,:));
    if mod(f,updateref)==0
        bead_im=IM_corr((bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),(bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),f);
        [bead_midpt_relxy,~,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
        bead_midpt_xy=round((bead_midpt_relxy-1)+(bead_midpt_xy-search_pad));
        
        ref(:,:,f)=IM_reduced(...
            (bead_midpt_xy(2)-im_pad):(bead_midpt_xy(2)+im_pad),...
            (bead_midpt_xy(1)-im_pad):(bead_midpt_xy(1)+im_pad),...
            f);
    else
        ref(:,:,f)=ref(:,:,f-1);
    end
end
loc=loc(:,[2,1]);
pos2=redfact*loc;


loc2=loc;
all_loc=[[loc1,(1:size(IM_reduced,3))'];[loc2,(1:size(IM_reduced,3))']];

%{
%Plot reference beads in all frames:
im3dscroll(ref,[],(0:size(ref,3)-1));colorbar;
colormap('gray');axis equal;caxis([0 300])
%}

%Plot reduced image with tracked position: less memory/faster scrolling
im3dscroll(IM_reduced,all_loc,(1:size(IM_reduced,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])

%Plot on original stack: more memory/slower scrolling
allpos=[[pos1,(1:size(IMstack,3))'];[pos2,(1:size(IMstack,3))']];
im3dscroll(IMstack,allpos,(0:size(IMstack,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])

%pos contains Bead2 position in the first column
temp=pos2;
temp(sorted,:)=pos2;
pos2=temp;

%Fourier filtering script:
%{
%Im is the image
V=Im(300,:)';
X=1:size(Im,2)';
%V=sin(pi*X)+sin(10*pi*X);

dL=1;
Fs=1./dL;

%{
%Interpolation onto a grid is not required here
Xq=min(X):dL:max(X);
if mod(numel(Xq),2)==1
    Xq=Xq(1:end-1);
end
Vq = interp1(X,V,Xq,'linear');
%}
Xq=X;
Vq=V;

L=numel(Xq);
F=fft(Vq);
P2 = abs(F/L);
P1 = P2(1:L/2+1);
f = Fs*(0:(L/2))/L;
fsym=[f(1:end-1),f(end-1:-1:1)];
F(fsym<0.005 | fsym>0.1)=0;

%plot(f,(P1).^2),xlim([0,5])


Vfilt=ifft(F);
figure(100),clf(100),
plot(Xq,Vq,'-b'),hold on
plot(Xq,(Vfilt),'-r'),hold on
[median(abs(Vfilt)),2.5*std((Vfilt))]
%}

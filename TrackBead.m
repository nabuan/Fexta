function [pos]=TrackBead(dirname)
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

%Parameters:
updateref=20;           %Steps after which reference image is recalculated
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
IMstack(:,:,1)=IM;
for f=2:numel(fname)
    IM=imread([dirname,fname(f).name]);
    IMstack(:,:,f)=IM;
end

IM_reduced=Reduce_Stack(IMstack,redfact,redfact,1,'Mean');
for f=1:size(IM_reduced,3)
    IM_reduced(:,:,f) = imgaussfilt(IM_reduced(:,:,f),2);
end
IM_corr=bsxfun(@minus,IM_reduced,mean(IM_reduced,3));

%Click on approximate center of the stage bead
figure(999);clf(999);
imshow(IM_corr(:,:,1),[]);
bead_midpt_xy=round(ginput(1));
bead_im=IM_corr((bead_midpt_xy(2)-search_pad):(bead_midpt_xy(2)+search_pad),(bead_midpt_xy(1)-search_pad):(bead_midpt_xy(1)+search_pad),1);
close(999)

%Re-adjust reference image to find bead
[tb_midpt_relxy,~,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
bead_midpt_xy=round((tb_midpt_relxy-1)+(bead_midpt_xy-search_pad));
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
        [tb_midpt_relxy,~,~]=imfindcircles(bead_im,round([0.5*im_pad im_pad]),'ObjectPolarity','dark');
        bead_midpt_xy=round((tb_midpt_relxy-1)+(bead_midpt_xy-search_pad));
        
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
%Plot reference beads in all frames:
im3dscroll(ref,[],(0:size(ref,3)-1));colorbar;
colormap('gray');axis equal;caxis([0 300])

%Plot reduced image with tracked position: less memory/faster scrolling
im3dscroll(IM_reduced,[loc,(1:size(IM_reduced,3))'],(1:size(IM_reduced,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])

%Plot on original stack: more memory/slower scrolling
im3dscroll(IMstack,[pos,(1:size(IMstack,3))'],(0:size(IMstack,3))');colorbar;
colormap('gray');axis equal;caxis([0 300])
%Contains Bead2 position in the first column
%--------------------------------------------------------------------------
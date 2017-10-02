%This script loads sample images and reduces image size using mean over
%local pixels. Reduction in size is not likely to affect precision of
%measurement(? test this claim empirically). Moreover, it will reduce
%variability and computational time in calculating normalized cross
%correlations used to track beads across images.

%Tests: 1. Try to detected Rois automatically
%Tests: 2. Check difference between high-res and low-res accuracy
%dirname='F:\Rohan\Dropbox\Lab\Other Projects\Fexta-DNA analysis\FEXTA\dat\Images\';
dirname='F:\Rohan\Dropbox\Lab\Other Projects\Fexta-DNA analysis\SameTypeBead\01DNA\dsDNA\';
%dirname='C:\Users\Naba\Dropbox\DATA\alphaSSB_ssDNA\28Sep2017\03DNA\t7Exo\t7Exo_1\';
fname=dir([dirname,'*.jpg']);
IM=imread([dirname,fname(1).name]);
IMall=zeros([size(IM),numel(fname)]);
IMall(:,:,1)=IM;
for f=2:numel(fname)
IM=imread([dirname,fname(f).name]);
IMall(:,:,f)=IM;
disp(f)
end

fact=4;
IMred=Reduce_Stack(IMall,fact,fact,1,'Mean');
% IMred=bsxfun(@minus,IMred,median(IMred,3));
% for f=1:size(IMred,3)
%     temp=IMred(:,:,f);
% end


BW=IMred;
%I1=abs(diff(IMred,1,1));I1=I1(:,2:end,:);
%I2=abs(diff(IMred,1,2));I2=I2(2:end,:,:);
%BW=(I1+I2);

%centers=cell(size(IMred,3),1);
r=[];
for f=1:size(IMred,3)
BW(:,:,f) = imgaussfilt(BW(:,:,f),2);
%[c, radii, metric] = imfindcircles(BW(:,:,f),[10 20]);
%r=[r;[c,f.*ones(size(c,1),1)]];
BW(:,:,f) = edge(BW(:,:,f),'sobel');
end
r=r(:,[2,1,3]);
rall=bsxfun(@times,r,[fact,fact,1]);
im3dscroll(BW,r,[1:size(BW,3)]);axis equal;colormap('gray');
im3dscroll(IMall,rall,[1:size(I1,3)]);axis equal;colormap('gray');
figure,plot(r(1:2:end,3),r(1:2:end,2),'o'),hold on;plot(r(2:2:end,3),r(2:2:end,2),'o'),hold on;
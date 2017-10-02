% This function reduces MatLab image stacks.
% Data format of Orig is preserved (uint8, uint16, etc).
% reduction_x, reduction_y, and reduction_z can be doubles, e.g. 2.1
% reduct_method is a string: 'Min', 'Max', 'Mean', 'Median', or 'SD'

function Im=Reduce_Stack(Orig,reduction_x,reduction_y,reduction_z,reduct_method)

sizeOrig=size(Orig);
if length(sizeOrig)==2
    sizeOrig(3)=1;
end
classOrig=class(Orig);

reduction_y=round(reduction_y);
if reduction_y<1
    reduction_y=1;
end
reduction_x=round(reduction_x);
if reduction_x<1
    reduction_x=1;
end
reduction_z=round(reduction_z);
if reduction_z<1
    reduction_z=1;
end

if reduction_z>sizeOrig(3)
    reduction_z=sizeOrig(3);
end
reduct_amount=[reduction_y,reduction_x,reduction_z];

sizeIm=fix(sizeOrig./reduct_amount);
[xx,yy,zz]=ind2sub(reduct_amount,1:prod(reduct_amount));
Im=zeros(sizeIm,classOrig);

if reduction_y==1 && reduction_x==1 && reduction_z==1
    Im=Orig;
elseif reduction_y>1 || reduction_x>1 || reduction_z>1
    for i=1:sizeIm(3)
        Orig_temp=Orig(1:reduct_amount(1)*sizeIm(1),1:reduct_amount(2)*sizeIm(2),1+(i-1)*reduct_amount(3):i*reduct_amount(3));
        
        Orig3=zeros([sizeIm(1:2),prod(reduct_amount)],classOrig);
        for ii=1:prod(reduct_amount)
            Orig3(:,:,ii)=Orig_temp(xx(ii):reduct_amount(1):size(Orig_temp,1),yy(ii):reduct_amount(2):size(Orig_temp,2),zz(ii));
        end
        
        if strcmp(reduct_method,'Max')
            Orig3 = max(Orig3,[],3);
        elseif strcmp(reduct_method,'Min')
            Orig3 = min(Orig3,[],3);
        elseif strcmp(reduct_method,'Mean')
            Orig3 = mean(double(Orig3),3);
        elseif strcmp(reduct_method,'Median')
            Orig3=sort(Orig3,3);
            if mod(size(Orig3,3),2)==0
                Orig3=mean(Orig3(:,:,size(Orig3,3)/2:size(Orig3,3)/2+1),3);
            else
                Orig3=Orig3(:,:,(size(Orig3,3)+1)/2);
            end
            %Orig3 = median(double(Orig3),3);
        elseif strcmp(reduct_method,'SD')
            Orig3 = std(double(Orig3),[],3);
        else
            error('Unknown method')
        end
        Im(:,:,i)=Orig3;
        display([num2str(i),'/',num2str(sizeOrig(3))])
    end
end
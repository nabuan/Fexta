function [Bead_1, Bead_2, I_R,I_ref] = Chose_Ref_Beads( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[filenm,pathnm]=uigetfile('*.jpg','Chose Reference Image');
I_ref=double(imread([pathnm,filenm]));

figure(1)
imshow(I_ref,[])
title('Chose the Bead on the trap')

rectL= getrect;

rectL=ceil(rectL);yL=rectL(1);xL=rectL(2);dyL=rectL(3);dxL=rectL(4);
Bead_1=I_ref(xL:xL+dxL,yL:yL+dyL);



title('Chose the Bead on the pipet')
rectS= getrect(figure(1));

rectS=ceil(rectS);yL=rectS(1);xL=rectS(2);dyL=rectS(3);dxL=rectS(4);
Bead_2=I_ref(xL:xL+dxL,yL:yL+dyL);

title('Chose the range of motion')

rectR= getrect;
rectR=ceil(rectR);yR=rectR(1);xR=rectR(2);dyR=rectR(3);dxR=rectR(4);
I_R=[xR,dxR,yR,dyR];
I_ref=I_ref(xR:xR+dxR,yR:yR+dyR);
close(figure(1));

end


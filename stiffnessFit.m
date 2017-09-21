
function [stiffness,StiffnessData]= findStiffnessFit(calib,PathName,stiffFile)
M=1:1:12;
L=length(M);
[Fst,stpath,FilterIndexStiff] = uigetfile({'*.*','All Files(*.*)'}, 'provide the Y file ');
b = strcat(stpath,Fst);
stiffFile=importdata(b);
n=find(all(stiffFile==0,2));
Stiffness=[];
PP=zeros(L,L+1);
for i=1:L
m=M(i);
ext(1,:)=stiffFile((1:n-1),2);
Fext(1,:)=stiffFile((1:n-1),3); 
Fext(2,:)=stiffFile((1:n-1),6);
Fext=0.17234*mean(Fext);


maxm=find(Fext==max(Fext));
minm=find(-Fext==max(-Fext));
Fext1=Fext(minm:maxm);
ext1=ext(minm:maxm);

p=polyfit(Fext1,ext1,m);
for(j=1:m+1)
PP(m,j)=1/p(j);
end
switch m
case 1
    
Fext1=Fext(minm+10:maxm-10);
ext1=ext(minm+10:maxm-10);

Y = @(x) p(1).*x +p(2);
Stiffness=[Stiffness;1/p(m)]
case 2
Y = @(x) p(1).*x.^2 + p(2).*x.^1 + p(3);
Stiffness=[Stiffness;1/p(m)]
case 3
Y = @(x) p(1).*x.^3 + p(2).*x.^2 + p(3).*x +p(4);
Stiffness=[Stiffness;1/p(m)]
case 4
Y = @(x) p(1).*x.^4 + p(2).*x.^3 + p(3).*x.^2 +p(4).*x.^1 +p(5);
Stiffness=[Stiffness;1/p(m)]
case 5
Y = @(x) p(1).*x.^5 + p(2).*x.^4 + p(3).*x.^3 +p(4).*x.^2 +p(5).*x.^1 + p(6);
Stiffness=[Stiffness;1/p(m)]
case 6
Y = @(x) p(1).*x.^6 + p(2).*x.^5 + p(3).*x.^4 + p(4).*x.^3 + p(5).*x.^2 +p(6).*x.^1 +p(7);
Stiffness=[Stiffness;1/p(m)]
case 7
Y = @(x) p(1).*x.^7 + p(2).*x.^6 + p(3).*x.^5 + p(4).*x.^4 + p(5).*x.^3 +p(6).*x.^2 +p(7).*x.^1 + p(8);
Stiffness=[Stiffness;1/p(m)]
case 8
Y = @(x) p(1).*x.^8 + p(2).*x.^7 + p(3).*x.^6 + p(4).*x.^5 + p(5).*x.^4 +p(6).*x.^3 +p(7).*x.^2 + p(8).*x + p(9);    
Stiffness=[Stiffness;1/p(m)]
case 9
Y = @(x) p(1).*x.^9 + p(2).*x.^8 + p(3).*x.^7 + p(4).*x.^6 + p(5).*x.^5 +p(6).*x.^4 +p(7).*x.^3 + p(8).*x.^2 + p(9).*x.^1+p(10);    
Stiffness=[Stiffness;1/p(m)]
case 10 
Y = @(x) p(1).*x.^10 + p(2).*x.^9 + p(3).*x.^8 + p(4).*x.^7 + p(5).*x.^6 +p(6).*x.^5 +p(7).*x.^4 + p(8).*x.^3 + p(9).*x.^2+p(10).*x.^1 + +p(11);    
Stiffness=[Stiffness;1/p(m)]
case 11 
Y = @(x) p(1).*x.^11 + p(2).*x.^10 + p(3).*x.^9 + p(4).*x.^8 + p(5).*x.^7 +p(6).*x.^6 +p(7).*x.^5 + p(8).*x.^4 + p(9).*x.^3+p(10).*x.^2 + +p(11).*x.^1++p(12);    
Stiffness=[Stiffness;1/p(m)]
case 12 
Y = @(x) p(1).*x.^12 + p(2).*x.^11 + p(3).*x.^10 + p(4).*x.^9 + p(5).*x.^8 +p(6).*x.^7 +p(7).*x.^6+ p(8).*x.^5 + p(9).*x.^4+ p(10).*x.^3 + +p(11).*x.^2++p(12).*x + p(13);    
Stiffness=[Stiffness;1/p(m)]
end

figure(i)
ext3=ext-Y(Fext); hold on
ext2=ext-Y(0) ;
ext4=Y(Fext)-Y(0);
plot(ext2,Fext,'mo','MarkerSize',2);hold on
plot(ext3,Fext,'r.')
plot(ext4,Fext,'k','linewidth',2)
title(m)

end


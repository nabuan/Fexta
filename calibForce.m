function [calbForce,minplat,maxplat,calib] = calibForce(Force,Ext,midpoint,tolerance) 

B = find(abs(diff(Force)) < tolerance);

C=find(Force(B)<0.3*max(Force));
B(C)=[];
B
a=B(1);
b=B(end);

c=find(diff(Ext)==0);

Force1=Force;
Force1(c)=[];
Ext(c)=[];

PP= polyfit(Ext(B),Force1(B),1);
X=Ext(a):0.01:Ext(b);

Y=PP(1).*X+PP(2);
%plot(X,Y,'r')
Y
Cal=midpoint /((Y(1)+Y(end))/2);
Force=Cal.*Force ;
calbForce=Force;
minplat=a;
maxplat=b;
calib=Cal;


end

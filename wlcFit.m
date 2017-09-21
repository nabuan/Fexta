function [FittedForce,FittedExt,df,dx,EWLC,FWLC] = wlcFit(Fext,ext)
ds(1)= 0.3419;
ds(2)= 0.0900;
ds(3)=  1.3615e+03;

WLC  = @(f) ds(1)*(1- 0.5*(ds(2)./f).^0.5 + (1/ds(3)).*f);

FWLC=0.1 : 0.01 :max(Fext) ;

EWLC= WLC(FWLC);
% axes(handles.DNA);
% plot(eWLC,FWLC,'r');hold on

str=find((Fext>=0.1)&(Fext<45));
c=str(1);
d=str(length(str));

xdata =@(Os) ext(c+1:d)+ Os(1);
Fstr=Fext(c+1:d);

xth= @(Os) WLC(Fstr+Os(2));
err = @(Os)trapz(Fstr+Os(2),(xdata(Os(1))- xth(Os)).^2);
Offset = fminsearch(err,[0.24 0.1]);
df=Offset(2);
dx=Offset(1);
FittedForce=Fext+df;
FittedExt=ext+dx;

end
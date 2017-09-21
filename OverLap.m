function [df,dx] = OverLap(Fext,ext,overlapVal)
ds(1)= 0.3419;
ds(2)= 0.0900;
ds(3)=  1.3615e+03;

WLC  = @(f) ds(1)*(1- 0.5*(ds(2)./f).^0.5 + (1/ds(3)).*f);

FWLC=0.1 : 0.01 :max(Fext) ;

EWLC= WLC(FWLC);

str=find((Fext>=0.1)&(Fext<45));
c=str(1);
d=str(length(str));

switch overlapVal
    
    case 1
 display('baseline & Xoffset')
xdata =@(Os) ext(c+1:d)+ Os(1);
Fstr=Fext(c+1:d);

xth= @(Os) WLC(Fstr+Os(2));
err = @(Os)trapz(Fstr+Os(2),(xdata(Os(1))- xth(Os)).^2);
Offset = fminsearch(err,[0.24 0.1]);
df=Offset(2);
dx=Offset(1);
    case 2
display('baseline only')  
xdata =ext(c+1:d);
Fstr=Fext(c+1:d);

xth= @(Os) WLC(Fstr+Os);
err = @(Os)trapz(Fstr+Os,(xdata- xth(Os)).^2);
Offset = fminsearch(err,1);
df=Offset;
dx=0;
    case 3
display('Xoffset only')  
xdata =@(Os) ext(c+1:d)+ Os;
Fstr=Fext(c+1:d);

xth= WLC(Fstr);
err = @(Os)trapz(Fstr,(xdata(Os)- xth).^2);
Offset = fminsearch(err,0.24);
df=0;
dx=Offset;
        
end        

end
function [Force,extension] = AnalayzeData(bp,Fdata,extdata,calibration,baseline,bL,p,xoffset) 



Fdata(1,:)=Fdata(1,:)-bL(1);
Fdata(2,:)=Fdata(2,:)-bL(2);
extdata=1000.*extdata./bp;


Fdata=calibration.*mean(Fdata) + baseline;
extdata=correctStiffness(p,extdata,Fdata,bp,0);

extdata=extdata+xoffset;


extension=extdata;

Force=Fdata;
end






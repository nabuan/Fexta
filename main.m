function [calibration,xoffset,blWLC,bL,p,EWLC,FWLC,minplat,maxplat] = main(bp,Fmp,feDNA,ext,stiffFile,tolerance,stiffOrder) 

midpoint=Fmp;
calibration=1;
blWLC=0;
xoffset=0;
ext=1000.*ext./bp;


% plot(ext,Fext,'r.','Parent',handles.DNA)

[feDNA(1,:),bs1]=Baseline(feDNA(1,:));
[feDNA(2,:),bs2]=Baseline(feDNA(2,:));
bL1=bs1;
bL2=bs2;

feDNA=mean(feDNA);

[fCal,minplat,maxplat,calib]=calibForce(feDNA,ext,midpoint,tolerance);

calibration=[calibration;calib];
p=findStiffnessFit(calib,stiffFile,stiffOrder)
ext=correctStiffness(p,ext,fCal,bp,0);



for(i=1:5)
[fCal,ext,bs,dx,EWLC,FWLC]=wlcFit(fCal,ext);
[fCal,minplat,maxplat,calib]=calibForce(fCal,ext,midpoint,tolerance);
calibration=[calibration;calib];
xoffset=[xoffset;dx];
blWLC=[blWLC;bs];
end

calibration=prod(calibration);
% stiffness=1/p(1);
xoffset=sum(xoffset);
blWLC=sum(blWLC);
bL=[bL1;bL2];
Force=fCal;
extension=ext;






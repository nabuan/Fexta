function calib = fitWLC(Dat,calib,Fmelt)

%Constants-----------------------------------------------------------------
kb=1.3806*10^-2;    %in pN*nm/K
T=298;              %in K
Pds=45;             %in nm
Bds=0.34;           %in nm/basePair
Sds=1361;           %in pN
Nbp=8100;           %User input - move this to parameter file later
%--------------------------------------------------------------------------

yoffset=calib.yoff;

%Conversion from um to nm.
Dat.xs=Dat.xs*10^3;
Dat.xt=Dat.xt*10^3;

xs=Dat.xs;
xt=Dat.xt;
xd=Dat.xd;

%ind=find(xd < (Fmelt-yoffset)./2 & xd > 9);
ind=find(xd < (Fmelt-yoffset)./1.5);

xs=xs(ind);
xt=xt(ind);
xd=xd(ind);

%Theoretical WLC: extension as a function of force
bds = @(F) Bds*(1-(1/2).* ...
    ((kb.*T)./(F.*Pds)).^0.5 + ...
    (F./Sds));

%figure,plot(X,Y)
%Ftemp = 0:0.1:60;
%figure,plot(bds(Ftemp),Ftemp)
% varmin = [c, xoffset]
argminf = @(varmin) ...
    sum(((xs-xt+varmin(2))./Nbp - ...
    Bds + ...
    (Bds./2)*((kb.*T)./((varmin(1).*(xd+yoffset-Fmelt)+Fmelt).*Pds)).^0.5 - ...
    (Bds./Sds).*(varmin(1).*(xd+yoffset-Fmelt) + Fmelt)).^2);

argminf2 = @(varmin) ...
    trapz(((varmin(1).*(xd+yoffset-Fmelt)+Fmelt)).^2,...
    ((xs-xt+varmin(2))./Nbp - ...
    Bds + ...
    (Bds./2)*((kb.*T)./((varmin(1).*(xd+yoffset-Fmelt)+Fmelt).*Pds)).^0.5 - ...
    (Bds./Sds).*(varmin(1).*(xd+yoffset-Fmelt) + Fmelt)).^2);

%bestvar=fminsearch(argminf,[0.1,400]);
bestvar=fminsearch(argminf2,[0.1,500]);

Fcalibrated=bestvar(1).*(Dat.xd+yoffset-Fmelt)+Fmelt;
bcalibrated=(Dat.xs-Dat.xt+bestvar(2))./Nbp;
Ftemp=[0.1:0.1:30];

h=gobjects(2);
figure(207);clf(207);
h(1)=plot(bcalibrated,Fcalibrated,'.','Color',[0.6 0.6 0.6],'MarkerSize',10);hold on
h(2)=plot(bds(Ftemp),Ftemp,'-','Color',[0.8 0 0]);
h(3)=plot(bcalibrated(ind),Fcalibrated(ind),'x','Color',[0 0 0],'MarkerSize',5);hold on
legend({'Calibrated data','WLC fit','Fitted points'},'Location','southeast');
grid on;axis square;drawnow;box on;
MakePretty(gca)

end
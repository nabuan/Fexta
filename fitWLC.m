function calib = fitWLC(Dat,calib,Fmelt)

%Constants-----------------------------------------------------------------
kb=1.3806*10^-2;    %in pN*nm/K
T=298;              %in K
Pds=45;             %in nm
Bds=0.34;           %in nm/basePair
Sds=1361;           %in pN
Nbp=48500;          %User input - move this to parameter file later
%--------------------------------------------------------------------------

yoffset=calib.yoff;

%Conversion from um to nm.
Dat.xs=Dat.xs*10^3;
Dat.xt=Dat.xt*10^3;

%Normalizing data using sigmoid fit parameters. Fmelt has to be changed as
%well, because we want to hold the y-position at the plateau at known
%value.
[f,ind] = fitSigmoid(Dat.xs-Dat.xt,Dat.xd); % -----> Estimate the sigmoid and choose datapoints to fit.
temp = (Dat.xd+yoffset-Fmelt)./f.A;
rough_y_scale = Fmelt./abs(min(temp(ind)));

%The following normalizes Dat.xd such that minimum is at 0 and the plateau
%is at Fmelt:
Dat.xd = ((Dat.xd+yoffset-Fmelt)./f.A)*rough_y_scale + Fmelt;
yoffset = 0;

Dat.xs=Dat.xs-f.B+0.35*Nbp;

xs=Dat.xs;
xt=Dat.xt;
xd=Dat.xd;

xs=xs(ind);
xt=xt(ind);
xd=xd(ind);

%Theoretical WLC: extension as a function of force
bds = @(F) Bds*(1-(1/2).* ...
    ((kb.*T)./(F.*Pds)).^0.5 + ...
    (F./Sds));

% argminf = @(varmin) ...
%     sum(...
%     ((xs-xt)./Nbp + varmin(2) -  ...
%     Bds + ...
%     (Bds./2)*((kb.*T)./((varmin(1).*(xd+yoffset-Fmelt)+Fmelt).*Pds)).^0.5 - ...
%     (Bds./Sds).*(varmin(1).*(xd+yoffset-Fmelt) + Fmelt)).^2);

argminf = @(varmin) ...
    sum(((xs-xt)./Nbp + varmin(2) - bds(varmin(1).*(xd+yoffset-Fmelt)+Fmelt)).^2);

%{
argminf2 = @(varmin) ...
    trapz(((varmin(1).*(xd+yoffset-Fmelt)+Fmelt)).^2,...
    ((xs-xt+varmin(2))./Nbp - ...
    Bds + ...
    (Bds./2)*((kb.*T)./((varmin(1).*(xd+yoffset-Fmelt)+Fmelt).*Pds)).^0.5 - ...
    (Bds./Sds).*(varmin(1).*(xd+yoffset-Fmelt) + Fmelt)).^2);
%}
    
opt = optimset('TolFun',10^-10,'TolX',10^-10);
[best_var,current_fval,exitflag,output]=fminsearch(argminf,[0.9,0],opt);
disp(best_var)

Fcalibrated=best_var(1).*(Dat.xd+yoffset-Fmelt)+Fmelt;
bcalibrated=(Dat.xs-Dat.xt)./Nbp+best_var(2);
Ftemp=(0.1:0.1:30);

h=gobjects(3,1);
figure(207);clf(207);
h(1)=plot(bcalibrated,Fcalibrated,'.','Color',[0.6 0.6 0.6],'MarkerSize',10);hold on
h(2)=plot(bds(Ftemp),Ftemp,'-','Color',[0.8 0 0]);
h(3)=plot(bcalibrated(ind),Fcalibrated(ind),'x','Color',[0 0 0],'MarkerSize',5);hold on
legend(h,{'Calibrated data','WLC fit','Fitted points'},'Location','southeast');
xlabel('Extension, nm/bp')
ylabel('Force, pN')
grid on;axis square;drawnow;box on;
MakePretty(gca)

end
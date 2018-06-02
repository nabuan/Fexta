function Yoffset = findYoffset(X,Y,F_melt)
%This function finds the approximate midpoint of the plateau region in the 
%data. The plateau corresponds to the DNA melting transition, and the
%midpoint of this transition can be used to calibrate the laser deflection.
%X=Dat.dsDNA.E.xs-Dat.dsDNA.E.xt;
%Y=Dat.dsDNA.E.xd;

X=X(:);Y=Y(:);
Xp=X(1:end-1);

dX=diff(X);
dY=diff(Y);
dYdX=dY./dX;

thr=quantile(dYdX,0.66);
indp=find(dYdX<thr);
ind=indp+1; % because numel(dYdX) = numel(X) - 1
ind=ind(Y(ind)>mean(Y));

X_melt=(X(ind(1))+X(ind(end)))./2;

%Linear regression to find Y_melt
ft = fittype('a*x+b');
[f,g] = fit(X(ind),Y(ind),ft,'StartPoint',[1 1]);
Y_melt=f(X_melt);
Yoffset=F_melt-Y_melt;

%{
disp(g)
h1=gobjects(3,1);
figure(206);clf(206)
h1(1)=plot(Xp,dYdX,'.-');hold on
h1(2)=plot([Xp(1),Xp(end)],[thr,thr],'-k');hold on
xlabel('Derivative, \delta{x}_{s}/\delta{x}_{d} (\mu{m})')
ylabel('Laser deflection x_{d} (A.U.)')
axis square,box on;grid on;drawnow;
MakePretty(gca);
%}

disp(g)
h=gobjects(4,1);
figure(205);clf(205)
h(1)=plot(X,Y,'.','Color',[0.7 0.7 0.7],'MarkerSize',10);hold on
h(2)=plot(X(ind),Y(ind),'.','Color',[0.2 0.2 0.2],'MarkerSize',10);hold on
h(3)=plot(X(ind),f(X(ind)),'-','Color',[0.8 0 0]);
h(4)=plot(X_melt,Y_melt,'^','Color',[0 0.4 0.8],'MarkerSize',10,'MarkerFaceColor',[0 0.4 0.8]);
legend(h,{'Data','Used to fit','Fit','Melting point reference'},'Location','southeast');
xlabel('Stage position, x_{s} (\mu{m})')
ylabel('Laser deflection x_{d} (A.U.)')
axis square,box on;grid on;drawnow;
MakePretty(gca);

end
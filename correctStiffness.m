function [xt] = correctStiffness(Expt,Stiff)
%This function returns position of the trapped bead for any deflection of
%the laser given the stiffness data.
%xt is the position of the trapped bead.
%Stiff is Dat.stiff.E;
%xd is average laser deflection in any other experiment
%Stiff.indmin and Stiff.indmax are necessary for the stiffness function to
%be invertible
xt = interp1(Stiff.xd(Stiff.indmin:Stiff.indmax),Stiff.xs(Stiff.indmin:Stiff.indmax),Expt.xd,'linear');
xt = xt - (- Stiff.fit.b./Stiff.fit.a);

%For debugging
figure(204),clf(204)
h(1)=plot(Expt.xs,Expt.xd,'.-','Color',[0.7 0.7 0.7],'MarkerSize',10);hold on
h(2)=plot(Expt.xs-xt,Expt.xd,'.-','Color',[0 0.8 0.4],'MarkerSize',10);hold on
legend(h,{'Stage position, x_{s}','Extention, x_{s} - x_{t}'})
xlabel('Position, (\mu{m})')
ylabel('Laser deflection x_{d} (A.U.)')
axis square,box on;grid on;drawnow;
MakePretty(gca);

end
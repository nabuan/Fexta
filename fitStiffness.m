function [f,indmin,indmax] = fitStiffness(xs,xd)
%Smooth data (using loess amoothing here)
%Find points where derivative changes value
%Restrict linear fitting to data within this range
%f is the fit function
%indmin and indmax are the derivative sign change locations

%Enforce sorted ordering of xs
[xs,sortind]=sort(xs);xd=xd(sortind);

%Detect points where the derivative changes sign
temp=sign(diff(smooth(xd,0.1,'loess')));
ind=find(diff(temp)~=0)+1;
[~,imin]=min(xd(ind));indmin=ind(imin);
[~,imax]=max(xd(ind));indmax=ind(imax);

%Linear region is calculated based on 50% of the xs range within these
%points (25% points removed from either side)
range_factor=0.25;
xs_range=(xs(indmax)-xs(indmin));

xsmin=xs(indmin)+range_factor*xs_range;
xsmax=xs(indmax)-range_factor*xs_range;
[~,indmin_lin]=min(abs(xs-xsmin));
[~,indmax_lin]=min(abs(xs-xsmax));

%Linear regression
ft = fittype('a*x+b');
[f,g] = fit(xs(indmin_lin:indmax_lin),xd(indmin_lin:indmax_lin),ft,'StartPoint',[1 1]);

%For debugging
disp(g)
h=gobjects(3,1);
figure(203);clf(203)
h(1)=plot(xs,xd,'.','Color',[0.7 0.7 0.7],'MarkerSize',10);hold on
h(2)=plot(xs(indmin_lin:indmax_lin),xd(indmin_lin:indmax_lin),'.','Color',[0.2 0.2 0.2],'MarkerSize',10);hold on
h(3)=plot(xs(indmin:indmax),f(xs((indmin:indmax))),'-','Color',[0.8 0 0]);
legend(h,{'All','Used to fit','Fit'})
xlabel('Stage position, x_{s} (\mu{m})')
ylabel('Laser deflection x_{d} (A.U.)')
axis square,box on;grid on;drawnow;
MakePretty(gca);
%}
end
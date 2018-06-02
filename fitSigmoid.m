function [fitind] = fitSigmoid(X,Y)

sig=fittype('A./(1+exp(-(x-B)/C))+D');

%Initial attempt
[sigfit,~]=fit(X,Y,sig,'StartPoint',[(max(Y)-min(Y))*0.75,mean(X),(max(X)-min(X))*0.20,min(Y)]);
fitind = find(X>(sigfit.B - 9*sigfit.C) & X<(sigfit.B + 9*sigfit.C));

%Refined fit
[sigfit,~]=fit(X(fitind),Y(fitind),sig,'StartPoint',[sigfit.A,sigfit.B,sigfit.C,min(Y(fitind))]);
fitind = find(X>(sigfit.B - 9*sigfit.C) & X<(sigfit.B));

figure(206);clf(206)
h1 = plot(X,Y,'x','Color',[0.4 0.4 0.4]);hold on;
h2 = plot(X,sigfit(X),'-r');
h3 = plot(X(fitind),Y(fitind),'o','Color',[0 0.5 0.8],'LineWidth',2);hold on;
grid on;box on; axis square;
ylabel('Force, A.U.')
xlabel('Extension, A.U.')
legend([h1,h2,h3],'Data','Sigmoid Fit','To fit with WLC')
MakePretty(gca);
end
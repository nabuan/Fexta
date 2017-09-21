
close all;
clear;
[DataFile,stpath,FilterIndexFile] = uigetfile({'*.*','All Files(*.*)'}, 'ForceFeedback data ');
DataF = strcat(stpath,DataFile);
data=importdata(DataF);

tData=data(:,1);
xData=data(:,2);
FData=data(:,3);

Force = mean(FData);
ErrOfForce =std(FData)/sqrt(length(Force));
tData=tData - min(tData);




Fit = @(pa,t) pa(1).*exp(-t./pa(2))+pa(3);
res = @(pa) sum(((Fit(pa,tData)-xData)).^2);
PA = fminsearch(res,[0.5,10,0.4]);
deltaX=PA(1)-PA(3);
str1= sprintf('x = %0.4f + %0.4fe^{-t/%0.4f} \n\ntau = %0.4f',PA(3),PA(1),PA(2),PA(2));

Fit2 = @(dpa,t) dpa(1).*(exp(-t./dpa(2)) + exp(-t./dpa(4)) )+dpa(3);
res2 = @(dpa) sum(((Fit2(dpa,tData)-xData)).^2);
DPA = fminsearch(res2,[0.5,10,0.4,1]);
deltaX=DPA(1)-DPA(3);

str2= sprintf('x = %0.4f + %0.4f ( e^{-t/%0.4f}+e^{-t/%0.4f} )\n\ntau1 = %0.4f \t\t tau2 = %0.4f',DPA(3),DPA(1),DPA(2),DPA(4),DPA(2),DPA(4));
str3=sprintf('@ Force = %0.4f +/- %0.4f pN', Force,ErrOfForce);
plot(tData,xData,'go');hold on
plot(tData,Fit(PA,tData),'b','LineWidth',2);
plot(tData,Fit2(DPA,tData),'r','LineWidth',2);

title(DataFile)
annotation('textbox',[0.35 0.75 0.50 0.075],'String',str3);
annotation('textbox',[0.35 0.60 0.50 0.15],'String',str1);
annotation('textbox',[0.35 0.45 0.50 0.15],'String',str2);
xlabel('time (s)')
ylabel('extension(nm/bp)')
% \n \delta x = %0.4f
% \n \delta = %0.4f

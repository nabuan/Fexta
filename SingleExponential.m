for i=1:10

[DataFile,stpath,FilterIndexFile] = uigetfile({'*.*','All Files(*.*)'}, 'ForceFeedback data ');
DataF = strcat(stpath,DataFile);
data=importdata(DataF);

tData=data(:,1);
xData=data(:,2);
FData=data(:,3);

Force = mean(FData);
ErrOfForce =std(FData)/sqrt(length(Force));
tData=tData - min(tData);


plot(tData,xData,'go');hold on

Fit = @(pa,t) pa(1).*exp(-t./pa(2))+pa(3);
res = @(pa) sum(((Fit(pa,tData)-xData)).^2);
PA = fminsearch(res,[0.5,10,0.4]);
deltaX=PA(1)-PA(3);
plot(tData,Fit(PA,tData),'b','LineWidth',2);
str= sprintf('x = %0.4f + %0.4fe^{-t/%0.4f} \n\nForce = %0.4f +/- %0.4f\n\ntau = %0.4f',PA(3),PA(1),PA(2),Force,ErrOfForce,PA(2));

title(DataFile)
annotation('textbox',[0.35 0.45 0.35 0.25],'String',str);
xlabel('time (s)')
ylabel('extension(nm/bp)')
% \n \delta x = %0.4f
% \n \delta = %0.4f
pause(10)
end
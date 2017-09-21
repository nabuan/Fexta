function stiffness= findStiffness(calib,stiffFile)
% cdold=pwd;
% cd(PathName);

n=find(all(stiffFile==0,2));

ext(1,:)=stiffFile((1:n-1),2);
Fext(1,:)=stiffFile((1:n-1),3); 
Fext(2,:)=stiffFile((1:n-1),6);
Fext=calib.*mean(Fext);


%plot(ext,Fext,'r.');hold on

maxm=find(Fext==max(Fext));
minm=find(-Fext==max(-Fext));
Fext1=Fext(minm:maxm);
ext1=ext(minm:maxm);
ind=find((Fext1>=-70)&(Fext1<70));
X=ext1(ind);
Y=Fext1(ind);
P = polyfit(X,Y,1) ;
stiffness=P(1);
XX=ext-Fext./stiffness;
StiffnessData(1,:)=XX;
StiffnessData(2,:)=Fext;


% % figure(3)
% str=sprintf('Stiffness=%f',stiffness);
% % plot(XX,Fext,'b.');
% title(str)


end
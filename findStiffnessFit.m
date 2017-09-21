
function p= findStiffnessFit(calibration,stiffFile,stiffOrder)


n=find(all(stiffFile==0,2));

ext(1,:)=stiffFile((1:n-1),2);
Fext(1,:)=stiffFile((1:n-1),3); 
Fext(2,:)=stiffFile((1:n-1),6);
Fext=calibration*mean(Fext);


maxm=find(Fext==max(Fext));
minm=find(-Fext==max(-Fext));
Fext1=Fext(minm:maxm);
ext1=ext(minm:maxm);

p=polyfit(Fext1,ext1,stiffOrder)







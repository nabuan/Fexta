
function SeperateForce(t,x,f,w,pathstr,DataName,OrgDir) 
tf=flip(t);
xf=flip(x);
ff=flip(f);

WindowSize=w;
r=filter(ones(1,WindowSize),1,t)./WindowSize;
rf=filter(ones(1,WindowSize),1,tf)./WindowSize;
r=r(WindowSize:end);
rf=rf(WindowSize:end);
u=filter(ones(1,WindowSize),1,x)./WindowSize;
uf=filter(ones(1,WindowSize),1,xf)./WindowSize;
u=u(WindowSize:end);
uf=uf(WindowSize:end);

v=filter(ones(1,WindowSize),1,f)./WindowSize;
vf=filter(ones(1,WindowSize),1,ff)./WindowSize;
v=v(WindowSize:end);
vf=vf(WindowSize:end);


figure(2)

% plot(t,x,'--');hold on
 plot(t,f,'--');hold on
%plot(r,v)

B = find(abs(diff(v)./diff(r))>0.01);
Bf = find(abs(diff(vf)./diff(rf))>0.01);


t2=t;
f2=f;
tf2=tf;
ff2=ff;
t2(B)=[];
f2(B)=[];
tf2(Bf)=[];
ff2(Bf)=[];

% plot(t2,f2,'ro')
% plot(tf2,ff2,'go')
C = find(abs(diff(f2))>1);
Cf = find(abs(diff(ff2))>1);

L=length(f);
n=length(C);
ind=zeros(1,n);
indf=zeros(1,n);
for i=1:n
   A=find(t==t2(C(i)+1));
   ind(i)=A(end);
   Af= find(t==tf2(Cf(i)+1));
    indf(i)=Af(end);
end
    
indf=flip(indf);

ind1=[1,ind];
indf1=[indf, L];
t_final=zeros(1,length(t),n);
f_final=zeros(1,length(t),n);
x_final=zeros(1,length(t),n);

for j=1:n+1
  
A=[t(ind1(j):indf1(j)),x(ind1(j):indf1(j)),f(ind1(j):indf1(j))]; 
assignin('base',sprintf('t_%d',j),A)
FileName=sprintf('%s_FF@%0.1fpN-anlyzd.txt',DataName,mean(A(:,3)));
cd(pathstr)
dlmwrite(FileName,A,'delimiter','\t','newline','pc');
[AX,H1,H2]=plotyy(A(:,1),A(:,3),A(:,1),A(:,2));hold on
V=get(H1,'color');
set(H1,'marker','o','markersize',3);
set(H2,'color',V,'marker','o','markersize',3);
set(AX(1),'Ylim',[0 max(f)],'YTick',[0:10:max(f)+5]);
set(AX(2),'Ylim',[min(x)-0.005 max(x)+0.005],'YTick',[min(x)-0.005:0.01:max(x)+0.005]);
 ylabel(AX(1),'Force (pN)')
    ylabel(AX(2),'Extension (nm/bp)')
    xlabel('time (s)')
end

% axes(AX(2))
% plot(t,x,'--');hold on

plot(t,f,'--');hold on
% set(AX(1),'Ylim',[0 max(f)],'YTick',[0:10:max(f)+5]);
% set(AX(2),'Ylim',[min(x)-0.005 max(x)+0.005],'YTick',[min(x)-0.005:0.01:max(x)+0.005]);


% plot(t(1:indf(1)),f(1:indf(1)),'yo')
% plot(t(ind(end):end),f(ind(end):end),'yo')
% % plot(t(ind(1):indf(2)),f((ind(1):indf(2))),'yo')
% % 
% % plot(t(ind(2):indf(3)),f((ind(2):indf(3))),'yo')
% 
% % plot(t2(C+1),f2(C+1),'m*')
% % plot(flip(tf2(Cf+1)),flip(ff2(Cf+1)),'k*')
%     m=n+1;
% while 1
% plot(t_1(:,1),t_1(:,2),'yo');hold on
% if m==1
%     break;
% end
% plot(t_2(:,1),t_2(:,2),'yo');hold on
% if m==2
%     break;
% end
% plot(t_3(:,1),t_3(:,2),'yo');hold on
% if m==3
%     break;
% end
% plot(t_4(:,1),t_4(:,2),'yo');hold on
% if m==4
%     break;
% end
% plot(t_5(:,1),t_5(:,2),'yo');hold on
% if m==5
%     break;
% end
% plot(t_6(:,1),t_6(:,2),'yo');hold on
% if m==6
%     break;
% end
% plot(t_7(:,1),t_7(:,2),'yo');hold on
% if m==7
%     break;
% end
% plot(t_8(:,1),t_8(:,2),'yo');hold on
% if m==8
%     break;
% end
% plot(t_9(:,1),t_9(:,2),'yo');hold on
% if m==9
%     break;
% end
% plot(t_10(:,1),t_10(:,2),'yo');hold on
% break;
% end
cd(OrgDir)
end
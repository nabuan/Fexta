Data=[1 1 1;0 0 0 ; 0 0 0; 1 1 1; 3 3 3;0 0 0 ; 5 5 5]

m=find(all(Data==0,2))

for(i=m(1):m(end))
if(all(Data(i+1,:)==0))
Data(i+1,:)=[];
break;
end
end

Data
m=find(all(Data==0,2))

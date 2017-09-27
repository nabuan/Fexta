function[StrRetnData,tData_FFE,tData_FFR] = SaveData(tData,xAData,FAData,imgIndData,ZeroRows)

type=length(ZeroRows);
m=ZeroRows;
k=length(tData);
tData_FFE=[];
tData_FFR=[];
switch(type)
   case 1
     Extension=1:1:m(1)-1;
     Return=m(1)+1:1:k;  
     Return=Return';
   case 2
      Extension=1:1:m(1)-1;
      FF_Extension=m(1)+1:1:m(2)-1;
      Return=m(2)+1:1:k;
      Return=Return';
    tData_FFE(:,1)=tData(FF_Extension);
    tData_FFE(:,2)=xAData(FF_Extension);
    tData_FFE(:,3)=FAData(FF_Extension);
    tData_FFE(:,4)=imgIndData(FF_Extension);
   case 3
      Extension=1:1:m(1)-1;
      Return1=m(1)+1:1:m(2)-1;
      FF_Return=m(2)+1:1:m(3)-1;
      Return3=m(3)+1:1:k;
      Return=[Return1';Return3'];
      
    tData_FFR(:,1)=tData(FF_Return);
    tData_FFR(:,2)=xAData(FF_Return);
    tData_FFR(:,3)=FAData(FF_Return);
    tData_FFR(:,4)=imgIndData(FF_Return);


  case 4
      Extension=1:1:m(1)-1;
      FF_Extension=m(1)+1:1:m(2)-1;
      Return1=m(2)+1:1:m(3)-1;
      FF_Return=m(3)+1:1:m(4)-1;
      Return3=m(4)+1:1:k;
      Return=[Return1';Return3']; 
      
    tData_FFE(:,1)=tData(FF_Extension);
    tData_FFE(:,2)=xAData(FF_Extension);
    tData_FFE(:,3)=FAData(FF_Extension);
    tData_FFE(:,4)=imgIndData(FF_Extension);

    tData_FFR(:,1)=tData(FF_Return);
    tData_FFR(:,2)=xAData(FF_Return);
    tData_FFR(:,3)=FAData(FF_Return);
    tData_FFR(:,4)=imgIndData(FF_Return);


end
teData=tData(Extension);
trData=tData(Return);
xeData=xAData(Extension);
xrData=xAData(Return);
FeData(1,:)=FAData(Extension);
FrData(1,:)=FAData(Return);
imgIndDataExt=imgIndData(Extension);
imgIndDataRtn=imgIndData(Return);

[xrData j]=sort(xrData);
FrData=FrData(j);
imgIndDataRtn=imgIndDataRtn(j);
StrRetnData=reshapeData([teData;xeData;FeData;imgIndDataExt],[trData;xrData;FrData;imgIndDataRtn]);

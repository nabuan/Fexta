function xData=correctStiffness(p,xrData,FrData,bp,STFILE)

xrData = (bp/1000).*xrData;

m=length(p)-1;
switch m
case 1
XofF_stiffness = @(x) p(1).*x +p(2);
case 2
XofF_stiffness = @(x) p(1).*x.^2 + p(2).*x.^1 + p(3);

case 3
XofF_stiffness = @(x) p(1).*x.^3 + p(2).*x.^2 + p(3).*x +p(4);

case 4
XofF_stiffness = @(x) p(1).*x.^4 + p(2).*x.^3 + p(3).*x.^2 +p(4).*x.^1 +p(5);

case 5
XofF_stiffness = @(x) p(1).*x.^5 + p(2).*x.^4 + p(3).*x.^3 +p(4).*x.^2 +p(5).*x.^1 + p(6);

case 6
XofF_stiffness = @(x) p(1).*x.^6 + p(2).*x.^5 + p(3).*x.^4 + p(4).*x.^3 + p(5).*x.^2 +p(6).*x.^1 +p(7);

case 7
XofF_stiffness = @(x) p(1).*x.^7 + p(2).*x.^6 + p(3).*x.^5 + p(4).*x.^4 + p(5).*x.^3 +p(6).*x.^2 +p(7).*x.^1 + p(8);

case 8
XofF_stiffness = @(x) p(1).*x.^8 + p(2).*x.^7 + p(3).*x.^6 + p(4).*x.^5 + p(5).*x.^4 +p(6).*x.^3 +p(7).*x.^2 + p(8).*x + p(9);    

case 9
XofF_stiffness = @(x) p(1).*x.^9 + p(2).*x.^8 + p(3).*x.^7 + p(4).*x.^6 + p(5).*x.^5 +p(6).*x.^4 +p(7).*x.^3 + p(8).*x.^2 + p(9).*x.^1+p(10);    

case 10 
XofF_stiffness = @(x) p(1).*x.^10 + p(2).*x.^9 + p(3).*x.^8 + p(4).*x.^7 + p(5).*x.^6 +p(6).*x.^5 +p(7).*x.^4 + p(8).*x.^3 + p(9).*x.^2+p(10).*x.^1 + +p(11);    

case 11 
XofF_stiffness = @(x) p(1).*x.^11 + p(2).*x.^10 + p(3).*x.^9 + p(4).*x.^8 + p(5).*x.^7 +p(6).*x.^6 +p(7).*x.^5 + p(8).*x.^4 + p(9).*x.^3+p(10).*x.^2 + +p(11).*x.^1++p(12);    

case 12 
XofF_stiffness = @(x) p(1).*x.^12 + p(2).*x.^11 + p(3).*x.^10 + p(4).*x.^9 + p(5).*x.^8 +p(6).*x.^7 +p(7).*x.^6+ p(8).*x.^5 + p(9).*x.^4+ p(10).*x.^3 + +p(11).*x.^2++p(12).*x + p(13);    

end

if(STFILE==1)
xData = xrData - XofF_stiffness(FrData);  
else
xData = xrData - XofF_stiffness(FrData)+XofF_stiffness(0);
end

xData = (1000/bp).*xData;

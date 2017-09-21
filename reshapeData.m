function [output] = reshapeData(A, B)

nA=length(A');
nB=length(B');
A=A';
B=B';
nA=length(A);
nB=length(B);
if(nA==nB)
    output=[A,B];
elseif (nA>nB)
    B(nB+1:nA,:)=NaN;
    output=[A,B];
elseif  (nA<nB)
    A(nA+1:nB,:)=NaN;
    output=[A,B];
end
end

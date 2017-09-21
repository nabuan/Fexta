x=0:1:10;
y1=x.^2;
y2=x;

[ax h1 h2]=plotyy(x,y1,x,y2);
delete(h2)
delete(h1)
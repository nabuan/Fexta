
%List of funtions for new version 

%Test for stiffness data:
%pth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/Stiffness.txt';
pth='C:\Users\Naba\Dropbox\Naba&Gala\FEXTA\dat\raw\Stiffness.txt';
vers=1;%default
typ='stiffness';
Dat=[];[Dat] = ReadFile(pth,typ,vers,Dat);

%pth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/DNA.txt';
pth='C:\Users\Naba\Dropbox\Naba&Gala\FEXTA\dat\raw\DNA.txt';
vers=1;%default
typ='dsDNA';
[Dat] = ReadFile(pth,typ,vers,Dat);

[Dat.stiff.E.fit,Dat.stiff.E.indmin,Dat.stiff.E.indmax] = fitStiffness(Dat.stiff.E.xs,Dat.stiff.E.xd);
[Dat.dsDNA.E.xt] = correctStiffness(Dat.dsDNA.E,Dat.stiff.E);
[F_melt] = ForceMidpoint(50);
[Dat.calibrate.yoff] = findYoffset(Dat.dsDNA.E.xs-Dat.dsDNA.E.xt,Dat.dsDNA.E.xd,F_melt);
[Dat.calibrate]=fitWLC(Dat.dsDNA.E,Dat.calibrate,F_melt);
%cftool with c./(1+(exp((x-e)./d)))+l./(1+(exp((x-m)./d)));
%This script reads stiffness and raw data, and returns the calibrated DNA curve

[~,val]=system('hostname');
if contains(val,'Fruity')
    %stiffnesspth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/Stiffness.txt';
    %rawpth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/DNA.txt';
    stiffnesspth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/stiffness,09,11,2017_5um.txt';
    rawpth='/Users/Fruity/Dropbox/Lab/Other Projects/Fexta-DNA analysis/FEXTA/dat/raw/dsDNA_48.5kbp_5umbead.txt';
else
    stiffnesspth='C:\Users\Nabuan\Dropbox\Naba&Gala\FEXTA\dat\raw\stiffness,09,11,2017_5um.txt';
    rawpth='C:\Users\Nabuan\Dropbox\Naba&Gala\FEXTA\dat\raw\dsDNA_48.5kbp_5umbead.txt';
end

vers=1;%default
typ='stiffness';
Dat=[];[Dat] = ReadFile(stiffnesspth,typ,vers,Dat);

vers=1;%default
typ='dsDNA';
[Dat] = ReadFile(rawpth,typ,vers,Dat);

[Dat.stiff.E.fit,Dat.stiff.E.indmin,Dat.stiff.E.indmax] = fitStiffness(Dat.stiff.E.xs,Dat.stiff.E.xd);
[Dat.dsDNA.E.xt] = correctStiffness(Dat.dsDNA.E,Dat.stiff.E);
[F_melt] = ForceMidpoint(50);
[Dat.calibrate.yoff] = findYoffset(Dat.dsDNA.E.xs-Dat.dsDNA.E.xt,Dat.dsDNA.E.xd,F_melt);
[Dat.calibrate]=fitWLC(Dat.dsDNA.E,Dat.calibrate,F_melt);
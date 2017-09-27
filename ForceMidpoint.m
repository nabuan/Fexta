function Fos = ForceMidpoint(SaltConc)

Na = [250,100,53.5,25,10,2.6]./1000;
fos= [65.5,62.6,61.0,58.8,55.3,51.5];
err= [0.8,0.7,0.7,1.2,1.2,1.2];
FOS = @(p,q) p(1).*log(q)+ p(2);
RES = @(p,q) sum(((FOS(p,Na)-fos)./err).^2);
P = fminsearch(RES,[2.5,60]);
Fos=FOS(P,SaltConc/1000);
end

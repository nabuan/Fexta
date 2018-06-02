function F_melt = ForceMidpoint(Na_Conc)
%Salt concentration Na_Conc input is in mM; 
%Calculations within this function are after conversion to M
%The hard-coded values are based on the following studies:
%Williams MC et al. (2001) ,Rouzina I et al. (2001) and Wenner RJ et al. (2002)

Na_Conc=Na_Conc/1000;
Na_ref = [250,100,53.5,25,10,2.6]./1000;
f_melt_theory= [65.5,62.6,61.0,58.8,55.3,51.5];
err= [0.8,0.7,0.7,1.2,1.2,1.2];
f_melt_fit = @(p,q) p(1).*log(q)+ p(2);
res = @(p,q) sum(((f_melt_fit(p,Na_ref)-f_melt_theory)./err).^2);
P = fminsearch(res,[2.5,60]);
F_melt=f_melt_fit(P,Na_Conc);
end

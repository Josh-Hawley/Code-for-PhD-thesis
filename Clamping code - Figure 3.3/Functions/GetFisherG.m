function [fisher_g,pval,idx]=GetFisherG(Pxx)
idx=find(Pxx==max(Pxx),1,'first');
fisher_g=Pxx(idx)/sum(Pxx);
N = length(Pxx);
    upper  = floor(1/fisher_g);
    for nn = 1:3
        I(nn) = (-1)^(nn-1)*nchoosek(N,nn)*(1-nn*fisher_g)^(N-1);
    end
pval = sum(I);

%GetFisherG is based on the matlab page and refs found here: 
%https://www.mathworks.com/help/signal/ug/significance-testing-for-periodic-component.html
function [fisher_g,pval,idx]=GetFisherG(Pxx)
    idx=find(Pxx==max(Pxx),1,'first');
    fisher_g=Pxx(idx)/sum(Pxx);
    % N = length(Pxx);
    % upper  = floor(1/fisher_g);
    %     for nn = 1:3
    %         I(nn) = (-1)^(nn-1)*nchoosek(N,nn)*(1-nn*fisher_g)^(N-1)
    %     end
    % pval = sum(I);
    N = length(Pxx);
    nn = 1:floor(1/fisher_g);

    I = (-1).^(nn-1).*exp(gammaln(N+1)-gammaln(nn+1)-gammaln(N-nn+1)).*(1-nn*fisher_g).^(N-1);
    
    pval = sum(I);
end
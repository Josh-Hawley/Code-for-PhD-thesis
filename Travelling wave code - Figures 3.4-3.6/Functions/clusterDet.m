%clusterDet.m Purpose is to remove synchronised cells from comparison
%matrix that are not adjacent to the diagonal values

function [Xout]=clusterDet(X)
N=size(X,1);
Xout=X;

for n=1:N
    
    [~,nf,~]=RunLength_M(X(n,n:end)); %nf is a vector that gives the running length of same value entries starting from n and goes forwards
    [~,nb,~]=RunLength_M(X(n,1:n));
    
    %Normal non-boundary cluster condition
    Xout(n,n+nf(1):end)=0;
    if n>1
        Xout(n,1:n-nb(end))=0;
    end
    
    %Boundary cluster condition
    if length(nb)==1 && X(n,end)==1 %If nb only has one element, this means that all vals to the left of the diagonal are = 1. Also if X on the right most inx has a value of one this means that the cluster wraps around through the periodic condition
%        nf(end)
%        Xout
        Xout(n,end-nf(end):end)=1; 
    end
    
end

for n=N:-1:2
    prior=sum(find(Xout(n-1,:)==1));
    current=sum(find(Xout(n,:)==1));
    if current==prior
        Xout(n,1:n-1)=0;
        Xout(n,n+1:end)=0;
    end
end
  

%KOP.m generates a Kuramoto order parameter for an array of phase angles 

function kop=KOPcomp(PH)

N=size(PH,1);
kop=zeros(N,N);

for i=1:N
    for j=1:N
        kop(i,j)=mean(abs(0.5*sum(exp(1i*[PH(i,:);PH(j,:)]))));   
    end
end



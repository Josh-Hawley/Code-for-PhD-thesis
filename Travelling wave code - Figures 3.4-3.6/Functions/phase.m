function [PH,PH_unwrap,COP]=phase(x,Y)

PH=zeros(numel(Y(:,1)),numel(x));
PH_unwrap=zeros(numel(Y(:,1)),numel(x));

for m=1:numel(Y(:,1))

    PH(m,:)=angle(hilbert(Y(m,:)));
    PH_unwrap(m,:)=unwrap(PH(m,:));    
       
end
COP=(1/numel(Y(:,1)))*sum(exp(1i*PH)); % Complex order parameter
end
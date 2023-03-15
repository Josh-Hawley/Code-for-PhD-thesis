function [Pnew]=improvedKymo(P,rows,cols,kymoCols)
Nt=size(P,2);
cells=rows*cols;
Pnew=P;

for ti=1:Nt
    
    Pgrid=reshape(P(:,ti),[rows,cols]);
    
    for n=2:2:rows
        
        if kymoCols==1
            Pgrid(n,:)=Pgrid(n,:);
        elseif kymoCols==2
            Pgrid(n,:)=(Pgrid(n,:) + [0, Pgrid(n,1:end-1)])/2; %take average of m-1 and m entries on n^th row
        end
        
    end
    
    Pnew(:,ti)=reshape(Pgrid,[cells, 1]);

end
end
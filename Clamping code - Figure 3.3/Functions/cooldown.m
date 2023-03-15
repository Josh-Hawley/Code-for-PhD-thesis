function [DIFF]=cooldown(DIFF,coolEl)

rows=size(DIFF,1);

for i=1:rows
    diff=DIFF(i,:);
    events=find(diff==1);
    diffNum=length(events);
    
    if diffNum>1
    
        for j=2:length(events)

            if events(j)-events(j-1)<coolEl
                DIFF(i,events(j))=0;
                
            end
        end
    end
end

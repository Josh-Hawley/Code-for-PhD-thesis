function [y,pos]=microClustDetect(vect,thresh)
% vect
p=polyfit([1:numel(vect)]',vect,3); 
f=polyval(p,1:numel(vect));
idx=find(f==max(f),1,'first');
if idx>1 && idx<numel(vect)
    % position at center
    kymo=vect;
    cent=kymo(idx);
    % fold around the max
    v1=kymo(1:idx);
    % flip around center
    v1=v1(end:-1:1);
    v2=kymo(idx:end);
    % fit a function to the average and estimate the band size
    n1=numel(v1);
    n2=numel(v2);
    if n1>n2
        vmean(1:n2)=(v1(1:n2)+v2(1:n2))/2;
        vmean(n2+1:n1)=v1(n2+1:n1);
    else
        vmean(1:n1)=(v1(1:n1)+v2(1:n1))/2;
        vmean(n1+1:n2)=v2(n1+1:n2);
    end
    lidx=find(vmean<thresh*max(vmean),1,'first');
    % size is double lidx since pattern is observed from center
    y=2*lidx-1;
    if ~isempty(lidx)& y>=2
        y=min(2*lidx-1,numel(vect));
        pos=idx;
    else
        y=NaN;
        pos=NaN;
    end
else % if center of pattern found at edge discard
    y=NaN;
    pos=NaN;
end
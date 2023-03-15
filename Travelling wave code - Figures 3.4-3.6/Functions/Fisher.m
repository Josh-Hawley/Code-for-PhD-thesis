function [FisherPer,Occurrence]=Fisher(P1,f, rows)

II=size(P1,2);

% KymLong=YDETREND2;
%     zs=zscore(KymLong); %z-scored data
%     [PXX2,F2]=periodogram(zs,[],[],1);         %Run periodogram and boostrap intervals
    PXX2=P1(2:end,:);
    F2=f(2:end);
    %Preallocation for storing stats
    fisherG=zeros(II,1);
    pval=zeros(II,1);
    peakLoc=zeros(II,1);
    peakHeight=zeros(II,1);
    FisherPer=zeros(II,1);

    for ii=1:II

        pxx=PXX2(:,ii); % periodogram
        
        [fisherG(ii),pval(ii),idx]=GetFisherG(pxx); % Find the peak and collect Fisher G-statistic
        peakLoc(ii)=idx;
        peakHeight(ii)=pxx(idx);
        if pval(ii)<0.05
            FisherPer(ii)=1/F2(idx);
        else
            FisherPer(ii)=nan;
        end

        if FisherPer(ii)==rows
            FisherPer(ii)=NaN;
        end


    end
    
    Occurrence=1-sum(isnan(FisherPer))/II;
    
end
function [Ydetrend, t, Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt)

    [Ydetrend, t, Ysmooth]=detrendSgolay(Y_raw, polyOrder, 0, T, frameLength);

    if frac==0
        startIDX=1;
    else
        startIDX=frac*Nt;
    end
    
    Ydetrend=Ydetrend(:,startIDX:end);
  
    Ysmooth=Ysmooth(:,startIDX:end);
    t=t(startIDX:end);
    
    %Fast Fourier transform (FT) of the detrended signal
    Y=fft(Ydetrend,[],2);
    L=length(Y);
    Fs=1/dt;
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y/L).^2;
    P1 = P2(:,1:L/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    
    %Calculating the average FT signal
    avgFourier=mean(P1,1);
%     if cells==1
%         avgFourier=P1;
%     else
%         sum(P1)/numel(P1(:,1),1);
%     end
    [~,i]=max(P1,[],2);
    ind_per=1./(60*f(i)); %Individual cell fourier dominant periods
    [~,I]=max(avgFourier); %Peak value index of the average FT (largest frequency contribution to signal)
    
    df=f(2)-f(1);
    tenP=0.1*f(I);
    tenP_idx=0.1*f(I)/df; %Ten percent in index value of peak frequency value in power spectrum, used in coherence calculation
%     figure(1)
%     clf
%     whos
%     subplot(311)
%     plot(Y_raw(1,:))
%     subplot(312)
%     plot(Ydetrend(1,:))
%     subplot(313)
%     plot(f,avgFourier)
    if I==1
        C1=0;
        f_C1=0;
        P_C1=avgFourier(I);
        
    elseif round(tenP_idx)==tenP_idx
        C1=trapz( f(I-tenP_idx:I+tenP_idx),avgFourier(I-tenP_idx:I+tenP_idx) );
        f_C1=f(I-tenP_idx:I+tenP_idx);
        P_C1=avgFourier(I-tenP_idx:I+tenP_idx);
        
    else
        [f_C1, P_C1]=accurateC1(f,avgFourier,I,tenP,tenP_idx); %Interpolates between fourier values to get more accurate coherence area  
        C1=trapz(f_C1,P_C1);
    end
    
    C2=trapz(f,avgFourier);
    coherence=C1/C2;
    
end

function [f_C1, P_C1]=accurateC1(f,P,I,tenP,tenP_idx)
    

    idx_im=I-floor(tenP_idx);
    idx_om=I-ceil(tenP_idx);
    idx_ip=I+floor(tenP_idx);
    idx_op=I+ceil(tenP_idx);
    
    f_im=f(idx_im); %Frequency value at inner, minor (Subscripts: i=inner, m=minor)
    f_ip=f(idx_ip); %Frequency value at inner, plus (Subscripts: i=inner, p=plus)
    f_om=f(idx_om);
    f_op=f(idx_op);
    
    P_im=P(idx_im);
    P_ip=P(idx_ip); 
    P_om=P(idx_om);
    P_op=P(idx_op);
    
    f_mm=f(I)-tenP;
    f_mp=f(I)+tenP; 
    
    smidge_m=abs((f_mm-f_im))/abs((f_om-f_im));
    smidge_p=abs((f_mp-f_ip))/abs((f_op-f_ip));
    P_mm=P_om*smidge_m + P_im*(1-smidge_m);
    P_mp=P_op*smidge_p + P_ip*(1-smidge_p);
    
    f_C1=[f_mm, f(idx_im:idx_ip), f_mp];
    P_C1=[P_mm,  P(idx_im:idx_ip), P_mp ];
    
end
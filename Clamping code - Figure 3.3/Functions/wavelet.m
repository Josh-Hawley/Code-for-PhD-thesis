%wavelet.m Calculates the wavelet transform of a real time signal across a 
%predefined frequency range.

%x  = signal to be analysed (time should be column-wise)
%w0 = lower bound frequency to analyse in Hz
%wf = upper bound frequency to analyse in Hz
%wc = central frequency for wavelet (wc<1 gives good time resolution but poor frequency resolution)
%Nw = number of frequency elements to use in wavelet analysis
%Ntime = number of time elements to use in wavelet analysis

% wc value is the central frequency, typical values range from 0.5 to 4: 
%                         0.5   1   2   3   4    
% Better time resolution <-------------------> Better frequency resolution

function [WT,W,T]=wavelet(x,Nw,w0,wf,wc,Ntime,t)
w0=w0; %lower frequency bound
wf=wf; %final frequency
dw=(wf-w0)/Nw;
ti=round(length(x)/Ntime); %Time indices to evaluate wavelet transform at

% W=(w0:dw:wf); %Vector of frequency values to be used
W=logspace(log10(w0), log10(wf), Nw);
T=t(1:ti:end);   %Vector of time values to be used

NW=length(W);
NT=length(T);
WT=zeros(NW,NT); %Wavelet transform

for i=1:NW %Loop through frequencies
    
    s=1/W(i);
    
    if ~mod(i,floor(NW*2/100)) %Percentage counter
        clc
        fprintf('Wavelet progress: %.0f%% \n', i/NW*100);
    end
    
    for j=1:NT %Loop through time
        
        ts=t-T(j); %ts = shifted time to slide wavelet along the signal 
        std=sqrt(2)*s; %standard deviation of the wavelet
        Nstd=2;
        
        if T(j)<Nstd*std || t(end)-T(j)<Nstd*std
            WT(i,j)=NaN;
        else
            phi=(1/pi^(1/4)) * (exp(2*pi*1i*wc.*ts./s) - exp(-2*pi*wc^2/2)) .*exp(-ts.^2/(2*s^2));
            WT(i,j)=trapz(t,phi.*x)/s; %Dividing by s normalises peak wavelet height due to wavelet area changing with scaling
    
        end
    
    end
end

W=W*wc;
W(all(isnan(WT),2))=[];
WT(all(isnan(WT),2),:)=[];

end


%         W(i)
%         figure(5)
%         clf
%         subplot(211)
%         plot(t,x); hold on
%         plot(t,real(phi),'w')
%         subplot(212)
%         plot(abs(WT(i,:)).^2)
%         drawnow
%         pause(0.01)
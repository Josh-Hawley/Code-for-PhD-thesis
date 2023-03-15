function [Ydetrend, t, Ysmooth]=detrendSgolay(Y_raw, polyOrder, frac, T, frameLength)
    cells=length(Y_raw(:,1));

    if cells<1
        error('Reduce Osc threshold')
    end


    Nt=length(Y_raw(1,:));

    if frac==0
        clip=1;
    else
        clip=ceil(frac*(Nt));
    end

    t=T(:,clip:end);

    %% Remove fraction of time vector and frequency vector

    Y_clip=Y_raw(:,clip:end);


    %% Detrend signal so that it lies around 0
    if mod(frameLength,2)==0
        frameLength=frameLength+1;
    end

    Ysmooth=sgolayfilt(Y_clip,polyOrder,frameLength,[],2);
    Ydetrend = Y_clip-Ysmooth; % Detrended y data

end
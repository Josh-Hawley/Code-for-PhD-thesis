function [mm]=MovingMean(A,wl,t_step)
    %A is the vector of protein values in time
    %wl is the window length which the mean should be taken
    
    if t_step<=wl
        mm=mean(A(:,1:t_step),2); %Moving mean value
    else
        mm=mean(A(:,t_step-wl:t_step),2);
    end
    
end
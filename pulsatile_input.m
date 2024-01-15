function [Ue, Ui] = pulsatile_input(multi,V_stimE,V_stimI,T_stim,x,duration,step,offduration)
Ue = zeros(1,duration/step);
Ui = zeros(1,duration/step);
duration1 = duration-offduration; 

% deals with any machine error in calcuation
epsilon = 0.001;

% ue input
counter = 0;
for i = 1:duration1/step
    counter = counter + step;
    
    % anodic phase
    if counter >= step && counter < T_stim+step-epsilon
        Ue(i) = -V_stimE/multi;
    end
    
    % cathodic phase
    if counter >= T_stim + step - epsilon && counter < (multi+1)*T_stim + step - epsilon
        Ue(i) = V_stimE;
    end
    
    if counter >= (multi+1)*T_stim + step - epsilon && counter < (multi+1+x)*T_stim  + step - epsilon
        Ue(i) = 0;
    end
    
    if counter >= (multi+1+x)*T_stim -  epsilon
        counter = 0;
        Ue(i) = 0;
    end
    
    
end

%ui input

counter = 0;
for i = 1:duration1/step
    counter = counter + step;
    
    % cathodic phase
    if counter >= step && counter < T_stim+step-epsilon
        Ui(i) = V_stimI;
    end
    
    % anodic phase
    if counter >= T_stim + step - epsilon && counter < (multi+1)*T_stim + step - epsilon
        Ui(i) = -V_stimI/multi;
    end
    
    % neutral phase
    if counter >= (multi+1)*T_stim + step - epsilon && counter < (multi+1+x)*T_stim  + step - epsilon
        Ui(i) = 0;
    end
    
    if counter >= (multi+1+x)*T_stim -  epsilon
       counter = 0;
       Ui(i) = 0;
    end
    
    
end
end
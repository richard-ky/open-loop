function[time,W_IE,RE1,v_Ems,v_Ims,S_EIms,S_IEms,S_EEms,S_IIms,X_EIms,X_IEms,X_EEms,X_IIms,Apostms,Aprems,W_IEms,spt_Ems,leftover_S_EI0,leftover_S_IE0,leftover_S_II0,leftover_S_EE0] = neuron_network(duration,Recordingduration,step,N_E,N_I,plast_on,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,J_EI,J_EE,J_IE,J_II,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,spt_E0)  

%clear -vars
 
   


    % -- Neuron Parameters --
    mew_e = 20.8; %constant current to each E neuron
    sigma_e = 1; % Gaussian noise standard deviation for E neuron
    mew_i = 18; % constant current to each I neuron
    sigma_i = 3; % Gaussian noise standard deviation for I neuron

    Ntot = N_E + N_I; %Total Number of Neurons
    
    

    % -- Synaptic Parameters --
    C_EI = 0.3*Ntot; % I-to-E connectivity
    C_EE = Ntot; % E-to-E connectivity
    C_IE = 0.3*Ntot; % E-to-I connectivity
    C_II = Ntot; % I-to-I connectivity
    tau_LTP = 20; %LTP time constant ms
    tau_LTD = 22; %LTD time constant ms

    %Initialization for Recording Data during 1s
    W_IE = zeros(Recordingduration/step,1); % E-to-I synaptic weights Storage
    spike_E = zeros(Recordingduration/step,N_E);
    spike_I = zeros(Recordingduration/step,N_I);
    time = zeros(Recordingduration/step,1);
    W_IEfs = zeros(Recordingduration/step,size(W_IE0,2));

    
   
    ref_E = zeros(1,N_E);
    ref_I = zeros(1,N_I);    
    tau_E_m = 10; %Excitatory neuron membrane time constan
    tau_I_m = 10; %Inhibitory neuron membrane time constant

    % -- Run --    
    
    sample_duration = 20;%Duration for each short simulation in ms
    Synchrony = zeros(duration/sample_duration,1);
    time_syn = zeros(duration/sample_duration,1);
    


    for ii = 1:duration/sample_duration
        % run parameters
        comp_time = (ii-1)*sample_duration;
       
    
        % indexes
        a = 1 + (ii>=2)*(ii-1)*sample_duration/step;
        b = ii*sample_duration/step;
        

        
    
       [timem, v_Em, v_Im, S_EIm, S_IEm, S_EEm, S_IIm, X_EIm, X_IEm, X_EEm, X_IIm, Apostm, Aprem, W_IEm, spike_Em, spike_Im, ref_Em, ref_Im, synchronym, spt_Em] = ode_neuron_model(plast_on,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,mew_e,sigma_e,mew_i,sigma_i,J_EI,J_EE,J_IE,J_II,C_EI,C_IE,C_EE,C_II,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,ref_E,ref_I,tau_E_m,tau_I_m,comp_time,spt_E0);
        
       % recorded variables

       %Spiking timing of E neurons
        spike_time_E(a:b,:) = spt_Em(1:sample_duration/step,:);
       if ii >= 951
       %time stamps of simulation
        c = a-190000;
        d = b-190000;
        %xpp = (c:d)
        time(c:d,:) = timem(1:sample_duration/step,:) + (ii-1)*sample_duration;
        
        %Mean E-to-I synaptic weight
        W_IE(c:d,1) = mean(W_IEm(1:sample_duration/step,:),2);
        
        %Spike event for E neruons
        spike_E(c:d,:) = spike_Em(1:sample_duration/step,:);


        %Spike event for I neruons
        spike_I(c:d,:) = spike_Im(1:sample_duration/step,:);
        
        %Synchrony measure for E population based on membrane potential
        Synchrony(ii) = synchronym;
        time_syn(ii) = sample_duration*(ii);

        

        %E-to-I synaptic weights
        %W_IEfs(c:d,:) = W_IEm(1:sample_duration/step,:);
       end

        % generate intial condition for next run
        sample_end = sample_duration/step ;
        vE0 = v_Em(sample_end,:);
        vI0 = v_Im(sample_end,:);
        S_EI0 = S_EIm(sample_end,:);
        S_IE0 = S_IEm(sample_end,:);
        S_II0 = S_IIm(sample_end,:);
        S_EE0 = S_EEm(sample_end,:);
        X_EI0 = X_EIm(sample_end,:);
        X_IE0 = X_IEm(sample_end,:);
        X_II0 = X_IIm(sample_end,:);
        X_EE0 = X_EEm(sample_end,:);
        Apost0 = Apostm(sample_end,:);
        Apre0 = Aprem(sample_end,:);
        W_IE0 = W_IEm(sample_end,:);
        left_sample_end = sample_end - 5/step;
        leftover_S_EI = S_EIm(left_sample_end:sample_end,:);
        leftover_S_IE = S_IEm(left_sample_end:sample_end,:);
        leftover_S_EE = S_EEm(left_sample_end:sample_end,:);
        leftover_S_II = S_IIm(left_sample_end:sample_end,:);
        spt_E0 = spt_Em(sample_end,:);
        
        %Output Variables
        v_Ems = vE0;
        v_Ims = vI0;
        S_EIms = S_EI0;
        S_IEms = S_IE0;
        S_EEms = S_EE0;
        S_IIms = S_II0;
        X_EIms = X_EI0;
        X_IEms = X_IE0;
        X_EEms = X_EE0;
        X_IIms = X_II0;
        Apostms = Apost0;
        Aprems = Apre0;
        W_IEms = W_IE0;
        spt_Ems = spt_E0;
        leftover_S_EI0 = leftover_S_EI;
        leftover_S_IE0 = leftover_S_IE;
        leftover_S_II0 = leftover_S_II;
        leftover_S_EE0 = leftover_S_EE;

    end


% calculated the Kuramoto Order Parameter
step = 0.1;
t = [0.1:step:duration];
[RE] = kuramoto_syn(spike_time_E,t,step,duration,N_E);
RE1 = RE(190001:200000,1);


function [RE] = kuramoto_syn(sptime,t,step,duration,N)
tspE = sptime;
phi = zeros(size(t,2),N);

% make vector of phases
for neuron = 1:N
        second_spike = 0;
    for i = 1:(duration)/step-1
        phi(i,neuron) = 2*pi*(t(i) - tspE(i,neuron));

        if tspE(i+1,neuron) ~= tspE(i,neuron)            
            if second_spike == 1
                delt = tspE(i+1,neuron) - tspE(i,neuron);
                a = floor(tspE(i,neuron)/step);
                b = floor(tspE(i+1,neuron)/step);
                unnorm = phi(a:b,neuron);
                phi(a:b,neuron) = unnorm./delt;
            end
            second_spike = 1;
        end
    end
end
RE = abs(mean((exp(1i.*phi)),2));
end
end
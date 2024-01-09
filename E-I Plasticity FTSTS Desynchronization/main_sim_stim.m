clear;
clc;



%Specify Average Synaptic Strength
J_IEN = (150:50:600)'; %E-to-I set this value 
nJIE = size(J_IEN,1);

%Specify Number of Neurons
N_E = 400; % # of E Neurons
N_I = 100; % # of I Neurons

%Specify Stimulation Parameters
T_stim = 1; %Pulse width of a single FTSTS pulse 0.5-2 - Interval of 0.5
x_neutral = 7; %Interpulse interval - Determine frequency
VstimE = 20; %E Neurons current amplitude %10 - 50 Interval of 5
VstimI = 20; %I Neurons current amplitude %10 - 50 Interval of 5
dxlE = 2; %Determine the spread of injected current to E neurons
dxlI = 0.5*dxlE; %Determine the spread of injected current to I neurons

multiple = 1;
interval_time = 4*4*(x_neutral*T_stim+T_stim+multiple*T_stim);

%Specify Simulation Parameters
duration = interval_time*150; %Duration for small batch
Stimoffduration = 1*1000; %Stimulation off duration for recording data
step = 0.1; %step size
Ns = 40; % Number of short runs
Nsf = 1;% Number of separate batch of runs for data collection

for trial = 1:nJIE
    J_IE = J_IEN(trial,1);
    %Load the initial steady-state data for initilaization
    FileName = ['SS_JIE',num2str(J_IE),'.mat'];
    
    FolderName = 'Staedy State';
    File = fullfile(FolderName, FileName);
    load(File)

    %Load the steady state converged synaptic weights
    FileName = ['synaptic_weight_data_', num2str(J_IE), '.mat'];
    FolderName = 'Staedy State';
    File = fullfile(FolderName, FileName);
    load(File)


%Finding total number of synapses in each connnectivity type
 num_synapses_IE = max(max(S_key_IE)); %E-to-I
 num_synapses_EI = max(max(S_key_EI)); %I-to-E
 num_synapses_EE = max(max(S_key_EE)); %E-to-E
 num_synapses_II = max(max(S_key_II)); %I-to-I

%Convert the initial synaptic weights to correct format
    for k = 1:N_E
        for j = 1:N_I
            if S_key_IE(k,j) ~= 0
                    index = S_key_IE(k,j);
                    W_IE0(1,index) = W_IEf0(k,j);
            end
        end
        for j = 1:N_E
            if S_key_EE(k,j) ~= 0
                index = S_key_EE(k,j);
                W_EE0(1,index) = W_EEf0(k,j);
            end 
        end
    end

   for k = 1:N_I
        for j = 1:N_E
            % synaptic input from I to E : _(EI)
            if S_key_EI(k,j) ~=0  
                index = S_key_EI(k,j);
                W_EI0(1,index) = W_EIf0(k,j); 
            end
        end
        for j = 1:N_I
            if S_key_II(j,k) ~= 0
                index = S_key_II(j,k);
                W_II0(1,index) =  W_IIf0(j,k);
            end
        end
   end

%Simulations

    disp(['Data Collection # ',' JIE ',num2str(J_IE),'AmpE',num2str(VstimE),'AmpI',num2str(VstimI),'PW',num2str(T_stim)])

        [REf,timef,W_IEA] = long_sim(Ns,Nsf,duration,Stimoffduration,step,N_E,N_I,T_stim,x_neutral,VstimE,VstimI,dxlE,dxlI,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,J_EI,J_EE,J_IE,J_II,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,spt_E0);

    save(['Stim_data_analysis_',' JIE ',num2str(J_IE),'AmpE',num2str(VstimE),'AmpI',num2str(VstimI),'PW',num2str(T_stim), '.mat'],'REf','timef','W_IEA','-v7.3')
end

minute = toc/60

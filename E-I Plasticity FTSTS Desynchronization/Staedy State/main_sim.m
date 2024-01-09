clear;
clc;


tic;
%Specify Average Synaptic Strength
J_IEN = (150:50:600)'; %E-to-I - set this value 
nJIE = size(J_IEN,1);
J_EI = 100; %I-to-E fixed
J_EE = 50; %E-to-E fixed
J_II = 50; %I-to-I fixed



%Specify Number of Neurons
N_E = 400; % # of E Neurons
N_I = 100; % # of I Neurons



%Specify Simulation Parameters
duration = 20*1000; %Duration for small batch
Recordingduration = 1*1000; %Duration for recording data
step = 0.1; %step size
Ns = 5; % Number of short runs - Ns*20 seconds of simulation
Nsf = 1;% Number of separate batch of runs for data collection - Nsf*Ns*20 seconds of simulation total

%Synaptic Data
load('Synaptic_Connection.mat');
for trial = 1:nJIE
    [W_IEf0,W_EIf0,W_IIf0,W_EEf0] = synptic_weight(N_E,N_I,S_key_EE,S_key_II,S_key_EI,S_key_IE);
    J_IE = J_IEN(trial,1);
    save(['synaptic_weight_data_', num2str(J_IE),'.mat'],'W_IEf0','W_EIf0','W_IIf0','W_EEf0','S_key_EE','S_key_II','S_key_EI','S_key_IE','J_IE','J_EI','J_EE','J_II')




    %Finding total number of synapses in each connnectivity type
    num_synapses_IE = max(max(S_key_IE)); %E-to-I
    num_synapses_EI = max(max(S_key_EI)); %I-to-E
    num_synapses_EE = max(max(S_key_EE)); %E-to-E
    num_synapses_II = max(max(S_key_II)); %I-to-I

    %  %Initialization
    vE0 = 14*ones(1,N_E);
    vI0 = 14*ones(1,N_I);
    S_EI0 = zeros(1,N_E);
    S_IE0 = zeros(1,N_I);
    S_EE0 = zeros(1,N_E);
    S_II0 = zeros(1,N_I);
    X_IE0 = zeros(1,N_I);
    X_EI0 = zeros(1,N_E);
    X_II0 = zeros(1,N_I);
    X_EE0 = zeros(1,N_E);
    spt_E0 = 0;

    Apost0 = zeros(1,num_synapses_IE);
    Apre0 = zeros(1,num_synapses_IE);

    leftover_S_EI = zeros(5/step + 1 ,N_E);
    leftover_S_IE = zeros(5/step + 1 ,N_I);
    leftover_S_EE = zeros(5/step + 1 ,N_E);
    leftover_S_II = zeros(5/step + 1 ,N_I);

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

    disp(['Data Collection # ',' JIE ',num2str(J_IE)])
   
    [REf,timef,W_IEA] = long_sim(Ns,duration,Recordingduration,step,N_E,N_I,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,J_EI,J_EE,J_IE,J_II,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,spt_E0);
    
    save(['Steady_State_data_analysis_', 'JIE',num2str(J_IE),'.mat'],'REf','timef','W_IEA','-v7.3')

save(['SS_JIE',num2str(J_IE),'.mat'],'vE0','vI0','S_EI0','S_IE0','S_EE0','S_II0','X_EI0','X_IE0','X_EE0','X_II0','Apost0','Apre0','W_IE0','leftover_S_EI','leftover_S_IE','leftover_S_EE','leftover_S_II','spt_E0')
end
minute = toc/60

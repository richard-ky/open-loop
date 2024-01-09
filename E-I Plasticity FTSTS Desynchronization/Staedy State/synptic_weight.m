function [W_IEf0,W_EIf0,W_IIf0,W_EEf0] = synptic_weight(N_E,N_I,S_key_EE,S_key_II,S_key_EI,S_key_IE)


%Specify Initial Synaptic Weight
Weight_0  =1;


    %Finding total number of synapses in each connnectivity type
    num_synapses_IE = max(max(S_key_IE)); %E-to-I
    num_synapses_EI = max(max(S_key_EI)); %I-to-E
    num_synapses_EE = max(max(S_key_EE)); %E-to-E
    num_synapses_II = max(max(S_key_II)); %I-to-I

    
    %Initial uniform random synaptic weights
    W_IE0 = 0.001+ Weight_0*rand(1,num_synapses_IE);   %E-to-I
    W_EI0 = 0.001+Weight_0*rand(1,num_synapses_EI);   %I-to-E
    W_EE0 = 0.001+Weight_0*rand(1,num_synapses_EE);   %E-to-I
    W_II0 = 0.001+Weight_0*rand(1,num_synapses_II);   %I-to-I
    
    %Storing Initial Synaptic Weight Matrix 
    for k = 1:N_E
        for j = 1:N_I
            if S_key_IE(k,j) ~= 0
                    index = S_key_IE(k,j);
                    W_IEf0(k,j) = W_IE0(1,index); %E-to-I - Dimension - N_E by N_I
            end
        end
        for j = 1:N_E
            if S_key_EE(k,j) ~= 0
                index = S_key_EE(k,j);
                W_EEf0(k,j) = W_EE0(1,index); %E-to-E - Dimension - N_E by N_E
            end 
        end
    end

    for k = 1:N_I
        for j = 1:N_E
            if S_key_EI(k,j) ~=0  
                index = S_key_EI(k,j);
                W_EIf0(k,j) = W_EI0(1,index); %I-to-E - Dimension - N_I by N_E
            end
        end
        for j = 1:N_I
            if S_key_II(j,k) ~= 0
                index = S_key_II(j,k);
                W_IIf0(j,k) =  W_II0(1,index); %I-to-I - Dimension - N_I by N_I
            end
        end
    end
end
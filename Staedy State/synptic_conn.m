function [S_key_EE,S_key_II,S_key_EI,S_key_IE] = synptic_conn()
N_E = 400;
N_I = 100;
% -- Make Random Synaptic Conncetions ---
    %Synaptic connectivity probability
    epsilon_E = 0.15; %I-to-E connectivety probability
    epsilon_I = 0.15; %E-to-I connectivety probability
    epsilon_EE = 0.15; %E-to-E connectivety probability
    epsilon_II = 0.15; %I-to-I connectivety probability
    
    %Varibales for storing synaptic connections information
    S_key_IE = zeros(N_E,N_I); %E-to-I
    S_key_EI = zeros(N_I,N_E); %I-to-E
    S_key_EE = zeros(N_E,N_E); %E-to-E
    S_key_II = zeros(N_I,N_I); %I-to-I

    %--- I to E ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_I
                syn_count = syn_count + 1;
                S_key_EI(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EI(pre_neuron,post_neuron) = 0;
            end

        end
    end

    % --- E to I ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_E
                syn_count = syn_count + 1;
                S_key_IE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_IE(pre_neuron,post_neuron) = 0;
            end

        end
    end

    % --- E to E ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_EE && pre_neuron ~= post_neuron
                syn_count = syn_count + 1;
                S_key_EE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EE(pre_neuron,post_neuron) = 0;
            end

        end
    end



    % --- I to I ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_II && post_neuron ~= pre_neuron
                syn_count = syn_count + 1;
                S_key_II(pre_neuron,post_neuron) = syn_count;
            else
                S_key_II(pre_neuron,post_neuron) = 0;
            end

        end
    end

    save('Synaptic_Connection.mat','S_key_II','S_key_EE','S_key_IE','S_key_EI')
end
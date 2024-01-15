
function[REf,timef,W_IEAf] = long_sim(Ns,duration,Recordingduration,step,N_E,N_I,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,J_EI,J_EE,J_IE,J_II,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,spt_E0)
   
    for ii = 1:Ns

                plast_on = 0; %STDP OFF
                

            if ii > 1
                %Initialization from the last run
                vE0 = v_Ems;
                vI0 = v_Ims;
                S_EI0 = S_EIms(end,:);
                S_IE0 = S_IEms(end,:);
                S_EE0 = S_EEms(end,:);
                S_II0 = S_IIms(end,:);
                X_EI0 = X_EIms;
                X_IE0 = X_IEms;
                X_EE0 = X_EEms;
                X_II0 = X_IIms;
                Apost0 = Apostms;
                Apre0 = Aprems;
                W_IE0 = W_IEms; 
                leftover_S_EI = leftover_S_EI0;
                leftover_S_IE = leftover_S_IE0;
                leftover_S_EE = leftover_S_EE0;
                leftover_S_II = leftover_S_II0;
                spt_E0 = spt_Ems;
            end

            %Call neuron_network.m
            [time,W_IEA,RE,v_Ems,v_Ims,S_EIms,S_IEms,S_EEms,S_IIms,X_EIms,X_IEms,X_EEms,X_IIms,Apostms,Aprems,W_IEms,spt_Ems,leftover_S_EI0,leftover_S_IE0,leftover_S_II0,leftover_S_EE0] = neuron_network(duration,Recordingduration,step,N_E,N_I,plast_on,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,W_II0,W_EE0,J_EI,J_EE,J_IE,J_II,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,spt_E0);  

            if ii == 1
                
                REf = RE; %Kuromoto Order Parameter
                
                timef = time; %time
                W_IEAf = W_IEA; % Mean E-to-I Synaptic Weights
            else
                
                W_IEAf = [W_IEAf;W_IEA];
                timef = [timef;time];
                REf = [REf;RE];
                
            end
            if ii == Ns
                vE0 = v_Ems;
                vI0 = v_Ims;
                S_EI0 = S_EIms(end,:);
                S_IE0 = S_IEms(end,:);
                S_EE0 = S_EEms(end,:);
                S_II0 = S_IIms(end,:);
                X_EI0 = X_EIms;
                X_IE0 = X_IEms;
                X_EE0 = X_EEms;
                X_II0 = X_IIms;
                Apost0 = Apostms;
                Apre0 = Aprems;
                W_IE0 = W_IEms; 
                leftover_S_EI = leftover_S_EI0;
                leftover_S_IE = leftover_S_IE0;
                leftover_S_EE = leftover_S_EE0;
                leftover_S_II = leftover_S_II0;
                spt_E0 = spt_Ems;
            end
            disp(ii) %Display each run
            save(['SS_Last_data_JIE',num2str(J_IE),'.mat'],'vE0','vI0','S_EI0','S_IE0','S_EE0','S_II0','X_EI0','X_IE0','X_EE0','X_II0','Apost0','Apre0','W_IE0','leftover_S_EI','leftover_S_IE','leftover_S_EE','leftover_S_II','spt_E0')
    end
end
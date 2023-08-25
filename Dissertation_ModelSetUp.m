%%  Fit Each Model
clear all; clc

% homedir='Y:\EEG_Data\JACKSON\RewP Model\RewP Model Summer 2020';  addpath(homedir); cd(homedir);
homedir='G:\UNM_Data\RewP Model\ModalityReward\';  addpath(homedir); addpath([homedir,'spm12\spm12\']); cd(homedir);

%% Model Parameters
% All caps used here are BrainVision 64-channel
% [2 36 21 60 11 45 15]; % Fz FCz Cz CPz Pz POz Oz

% ##### Soft Coded Stuff ##### %

tx2disp=-500:2:1000;

% Set Plausible Param Ranges For All Models

% ########
model = {'Full','Half','Gauss','Wavelet','TwoHalf'};% Can be any combination of {'Full','Half','Gauss','Wavelet','Gamma','TwoHalf'};
method = 'Grad'; % Can be 'Grid' (Grid search) or 'Grad' (Gradient Descent)
RMRUNS = 200; % initial sample
options=optimset('display','off','LargeScale','off'); % rmsearch options

% Params for models
FREQ_VECTOR=            0.5:.1:4;
SCALING_VECTOR=         -20:.1:20;
TIMESHIFT_VECTOR=       0:1:200;
TAPER_VECTOR=           10:1:100;
SHAPE_VECTOR=           5:5:40;
RATE_VECTOR=            0.2:0.05:1;
TwoModel_SCALING_VECTOR = 0:0.1:20;
TwoModel_TIMESHIFT_VECTOR = 1:1:200;
% ########

% ############################ %

%% Run Models

for UseDataset = 2%1:6
    if UseDataset == 1
        % R01 Door
        Study = 'E1_SightTrain_R01Door\'; cd([homedir,Study]);
        %load('R01DoorVolt_30-Sep-2021.mat');
        load('R01DepDoorVolt_15-Aug-2022.mat');
        srate = 500;
        DataTitle = 'R01_Door_';
        CondiTitle = {'Sight'};
        Iterations = 1;
        LossCondi = 3; % Specific to this Doors Task. 2 = Neutral; 3 = Loss;
        Electrode = 36;
        N = size(DoorVolt.ERPs_FB,1);
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,LossCondi) );
        Training_Sight.Rew_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,1) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
        
    elseif UseDataset == 2
        % Sight Or Sound
        Study = 'E2_SightOrSound\'; cd([homedir,Study]);
        load('SoS_Data.mat','MEGA_ERPs');
        srate = 500;
        DataTitle = 'SightOrSound_';
        CondiTitle = {'Sight','Sound'};
        Iterations = 2;
        LossCondi = [4 1];
        N = size(MEGA_ERPs.FBbyCUE,1);
        Electrode = 36;
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi(1)) );
        Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
        
        % Individual differences
        Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi(2)) );
        Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3) );
        
        % Everything Else
        Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
        Training_Sound.time1 = 100;
        Training_Sound.time2 = 500;
        Training_Sound.UseBaseline = 0;
        % Training_Sound.BDI = DoorVolt.BDI';
        % Training_Sound.BIS = DoorVolt.BIS';
        % Training_Sound.BAS_D = DoorVolt.BAS_D';
        % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sound.age = DoorVolt.Age';
        % Training_Sound.sex = DoorVolt.sex';
        
    elseif UseDataset == 3
        
        % Sight And Sound
        Study = 'E3_SightAndSound\'; cd([homedir,Study]);
        load('SnS_Data.mat');
        srate = 500;
        DataTitle = 'SightAndSound_';
        CondiTitle = {'Sight','Sound','Combo'};
        Iterations = 3;
        LossCondi = 1;
        N = size(MEGA_ERPs.FBbyCUE,1);
        Electrode = 36;
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
        Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
        
        % Individual differences
        Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
        Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,4) );
        
        % Everything Else
        Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
        Training_Sound.time1 = 100;
        Training_Sound.time2 = 500;
        Training_Sound.UseBaseline = 0;
        % Training_Sound.BDI = DoorVolt.BDI';
        % Training_Sound.BIS = DoorVolt.BIS';
        % Training_Sound.BAS_D = DoorVolt.BAS_D';
        % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sound.age = DoorVolt.Age';
        % Training_Sound.sex = DoorVolt.sex';
        
        % Individual differences
        Training_Combo.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
        Training_Combo.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2) );
        
        % Everything Else
        Training_Combo.RewP_Diff = Training_Combo.Rew_ERP - Training_Combo.Pun_ERP;
        Training_Combo.time1 = 100;
        Training_Combo.time2 = 500;
        Training_Combo.UseBaseline = 0;
        % Training_Combo.BDI = DoorVolt.BDI';
        % Training_Combo.BIS = DoorVolt.BIS';
        % Training_Combo.BAS_D = DoorVolt.BAS_D';
        % Training_Combo.BAS_FS = DoorVolt.BAS_FS';
        % Training_Combo.BAS_RR = DoorVolt.BAS_RR';
        % Training_Combo.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Combo.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Combo.age = DoorVolt.Age';
        % Training_Combo.sex = DoorVolt.sex';
        
    elseif UseDataset==4
        
        % Gremlin
        Study = 'E4_SoundDelay_Gremlin\'; cd([homedir,Study]);
        load('V_Gremlin_16-Aug-2022.mat');
        srate = 500;
        DataTitle = 'V_Gremlin_16-Aug-2022';
        CondiTitle = {'100 ms','300 ms','500 ms'};
        Iterations = 3;
        LossCondi = 1;
        N = size(MEGA_ERPs.FBbyCUE,1);
        Electrode = 21;
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,1) );
        Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,1) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 600;
        Training_Sight.UseBaseline = 100;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
        
        % Individual differences
        Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,2) );
        Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,2) );
        
        % Everything Else
        Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
        Training_Sound.time1 = 300;
        Training_Sound.time2 = 800;
        Training_Sound.UseBaseline = 300;
        % Training_Sound.BDI = DoorVolt.BDI';
        % Training_Sound.BIS = DoorVolt.BIS';
        % Training_Sound.BAS_D = DoorVolt.BAS_D';
        % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sound.age = DoorVolt.Age';
        % Training_Sound.sex = DoorVolt.sex';
        
        % Individual differences
        Training_Combo.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,3) );
        Training_Combo.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,3) );
        
        % Everything Else
        Training_Combo.RewP_Diff = Training_Combo.Rew_ERP - Training_Combo.Pun_ERP;
        Training_Combo.time1 = 500;
        Training_Combo.time2 = 1000;
        Training_Combo.UseBaseline = 500;
        % Training_Sound.BDI = DoorVolt.BDI';
        % Training_Sound.BIS = DoorVolt.BIS';
        % Training_Sound.BAS_D = DoorVolt.BAS_D';
        % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sound.age = DoorVolt.Age';
        % Training_Sound.sex = DoorVolt.sex';
        
        elseif UseDataset == 5
        % R01 Depressed Door
        Study = 'T1_SightDep_R01Door\'; cd([homedir,Study]);
        %load('R01DoorVolt_30-Sep-2021.mat');
        load('R01DepDoorVolt_14-Mar-2023.mat');
        srate = 500;
        DataTitle = 'R01_Door_';
        CondiTitle = {'Sight'};
        Iterations = 1;
        LossCondi = 3; % Specific to this Doors Task. 2 = Neutral; 3 = Loss;
        Electrode = 36;
        N = size(DoorVolt.ERPs_FB,1);
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,LossCondi) );
        Training_Sight.Rew_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,1) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
        
    elseif UseDataset == 6
        
        % Cavanagh et al., 2019
        Study = 'T2_RL_Dep\'; cd([homedir,Study]);
        %load('R01DoorVolt_30-Sep-2021.mat');
        load('PS_RewSens.mat');
        srate = 500;
        DataTitle = 'PST_Depression_Cav2019_PEs_';
        CondiTitle = {'Lo PE','Hi PE'};
        Iterations = 2;
        LossCondi = 1; 
        WinCondi=2;
        Electrode = 36;
        DepressedParticipants=DEP;
        N = size(MEGA_ERP,1);
        
        % Pull info for model
        
        % Individual differences (low)
        Training_Sight.Pun_ERP = squeeze( MEGA_ERP_PEs(:,Electrode,1:751,LossCondi,1) );
        Training_Sight.Rew_ERP = squeeze( MEGA_ERP_PEs(:,Electrode,1:751,WinCondi,1) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        
        % Individual differences (hi)
        Training_Sound.Pun_ERP = squeeze( MEGA_ERP_PEs(:,Electrode,1:751,LossCondi,2) );
        Training_Sound.Rew_ERP = squeeze( MEGA_ERP_PEs(:,Electrode,1:751,WinCondi,2) );
        
        % Everything Else
        Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
        Training_Sound.time1 = 100;
        Training_Sound.time2 = 500;
        Training_Sound.UseBaseline = 0;
        
    end
    % Outputs must be: [Pun_ERP    Rew_ERP    RewP_Diff    tx2disp   srate   DataTitle]
    % Where each ERP is a 1D vector, tx2disp is the samples-to-time conversion, srate is sampling rate, DataTitle is for saving the output
    
    
    
    %% Fit_RewP to all dem models
    for Iteri = 1:Iterations
        
        if Iteri == 1
            Use = Training_Sight;
        elseif Iteri == 2
            Use = Training_Sound;
        elseif Iteri == 3
            Use = Training_Combo;
        end
        
        for Modeli = 1:size(model,2)
            for n = 1:N
                
                
                if strcmp(model(Modeli),'Full') && ~exist([homedir,Study,DataTitle,'FullCosine_ModelFits_',method,'.mat'])
                    % Full Cosine model
                    
                    if all(method == 'Grid')
                        tic
                        % Standard Grid Search
                        for FREQi=1:length(FREQ_VECTOR)
                            for SCALINGi=1:length(SCALING_VECTOR)
                                for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                                    params=[FREQ_VECTOR(FREQi),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR(TIMESHIFTi)];
                                    %                         Model_AIC_Full(FREQi,SCALINGi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                    Model_AIC_Full(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                    
                                    clear params;
                                end % time end
                            end % Scaling end
                        end % Frequency end
                        Fulltoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Full);
                        n_data = dims(1).*dims(2).*dims(3);
                        ReshapeDaMatrix = reshape(Model_AIC_Full,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D]=ind2sub(size(Model_AIC_Full),find(Model_AIC_Full==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_FULL(n,Iteri).AIC=lowest_AIC;
                        GRID_FULL(n,Iteri).FREQ=FREQ_VECTOR(B);
                        GRID_FULL(n,Iteri).SCALING=SCALING_VECTOR(C);
                        GRID_FULL(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR(D);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'FullCosine_ModelFits_',method],'GRID_FULL');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(FREQ_VECTOR),median(SCALING_VECTOR),median(TIMESHIFT_VECTOR)]';
                        lower_limits = [min(FREQ_VECTOR),min(SCALING_VECTOR),min(TIMESHIFT_VECTOR)]';
                        upper_limits = [max(FREQ_VECTOR),max(SCALING_VECTOR),max(TIMESHIFT_VECTOR)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_FULL(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest SSE (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_FULL(n,Iteri).AIC=AIC_out;
                        GRAD_FULL(n,Iteri).FREQ=params_out(1);
                        GRAD_FULL(n,Iteri).SCALING=params_out(2);
                        GRAD_FULL(n,Iteri).TIMESHIFT=params_out(3);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'FullCosine_ModelFits_',method],'GRAD_FULL');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                elseif strcmp(model(Modeli),'Half') && ~exist([homedir,Study,DataTitle,'HalfCosine_ModelFits_',method,'.mat'])
                    
                    if all(method == 'Grid')
                        tic
                        % Half Cosine model
                        for FREQi=1:length(FREQ_VECTOR)
                            for SCALINGi=1:length(SCALING_VECTOR)
                                for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                                    
                                    params=[FREQ_VECTOR(FREQi),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR(TIMESHIFTi)];
                                    % Model_AIC_Half(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                    Model_AIC_Half(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                    
                                    clear params;
                                end % Timeshift end
                            end % Scaling end
                        end % Frequency end
                        Halftoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Half);
                        n_data = dims(1).*dims(2).*dims(3);
                        ReshapeDaMatrix = reshape(Model_AIC_Half,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D]=ind2sub(size(Model_AIC_Half),find(Model_AIC_Half==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_HALF(n,Iteri).AIC=lowest_AIC;
                        GRID_HALF(n,Iteri).FREQ=FREQ_VECTOR(B);
                        GRID_HALF(n,Iteri).SCALING=SCALING_VECTOR(C);
                        GRID_HALF(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR(D);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'HalfCosine_ModelFits_',method],'GRID_HALF');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(FREQ_VECTOR),median(SCALING_VECTOR),median(TIMESHIFT_VECTOR)]';
                        lower_limits = [min(FREQ_VECTOR),min(SCALING_VECTOR),min(TIMESHIFT_VECTOR)]';
                        upper_limits = [max(FREQ_VECTOR),max(SCALING_VECTOR),max(TIMESHIFT_VECTOR)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_HALF(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest SSE (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_HALF(n,Iteri).AIC=AIC_out;
                        GRAD_HALF(n,Iteri).FREQ=params_out(1);
                        GRAD_HALF(n,Iteri).SCALING=params_out(2);
                        GRAD_HALF(n,Iteri).TIMESHIFT=params_out(3);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'HalfCosine_ModelFits_',method],'GRAD_HALF');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                elseif strcmp(model(Modeli),'Gauss') && ~exist([homedir,Study,DataTitle,'Gaussian_ModelFits_',method,'.mat'])
                    
                    TIMESHIFT_VECTOR_Gauss=       0:1:451;
                    
                    if all(method == 'Grid')
                        tic
                        % Gaussian model
                        for FREQi=1:length(FREQ_VECTOR)
                            for SCALINGi=1:length(SCALING_VECTOR)
                                for TIMESHIFTi=1:length(TIMESHIFT_VECTOR_Gauss)
                                    
                                    params=[FREQ_VECTOR(FREQi),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR_Gauss(TIMESHIFTi)];
                                    % Model_AIC_Gauss(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                    Model_AIC_Gauss(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                    
                                    clear params;
                                end % Timeshift end
                            end % Scaling end
                        end % Frequency end
                        Gausstoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Gauss);
                        n_data = dims(1).*dims(2).*dims(3);
                        ReshapeDaMatrix = reshape(Model_AIC_Gauss,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D]=ind2sub(size(Model_AIC_Gauss),find(Model_AIC_Gauss==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_GAUS(n,Iteri).AIC=lowest_AIC;
                        GRID_GAUS(n,Iteri).FREQ=FREQ_VECTOR(B);
                        GRID_GAUS(n,Iteri).SCALING=SCALING_VECTOR(C);
                        GRID_GAUS(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR_Gauss(D);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'Gaussian_ModelFits_',method],'GRID_GAUS');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(FREQ_VECTOR),median(SCALING_VECTOR),median(TIMESHIFT_VECTOR_Gauss)]';
                        lower_limits = [min(FREQ_VECTOR),min(SCALING_VECTOR),min(TIMESHIFT_VECTOR_Gauss)]';
                        upper_limits = [max(FREQ_VECTOR),max(SCALING_VECTOR),max(TIMESHIFT_VECTOR_Gauss)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_GAUS(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest SSE (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_GAUS(n,Iteri).AIC=AIC_out;
                        GRAD_GAUS(n,Iteri).FREQ=params_out(1);
                        GRAD_GAUS(n,Iteri).SCALING=params_out(2);
                        GRAD_GAUS(n,Iteri).TIMESHIFT=params_out(3);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'Gaussian_ModelFits_',method],'GRAD_GAUS');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                elseif strcmp(model(Modeli),'Wavelet') && ~exist([homedir,Study,DataTitle,'GaussianTaper_ModelFits_',method,'.mat'])
                    
                    TIMESHIFT_VECTOR_Wave = 0:1:451;
                    
                    if all(method == 'Grid')
                        tic
                        % Wavelet model
                        for FREQi=1:length(FREQ_VECTOR)
                            for SCALINGi=1:length(SCALING_VECTOR)
                                for TAPERSHIFTi=1:length(TAPER_VECTOR)
                                    for TIMESHIFTi=1:length(TIMESHIFT_VECTOR_Wave)
                                        
                                        params=[FREQ_VECTOR(FREQi),SCALING_VECTOR(SCALINGi),TAPER_VECTOR(TAPERSHIFTi),TIMESHIFT_VECTOR_Wave(TIMESHIFTi)];
                                        %                             Model_AIC_Wavelet(FREQi,SCALINGi,TAPERSHIFTi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                        Model_AIC_Wavelet(FREQi,SCALINGi,TAPERSHIFTi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                        
                                        clear params;
                                    end % Timeshift end
                                end % Tapershift end
                            end % Scaling end
                        end % Frequency end
                        Wavetoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Wavelet);
                        n_data = dims(1).*dims(2).*dims(3).*dims(4);
                        ReshapeDaMatrix = reshape(Model_AIC_Wavelet,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D,E]=ind2sub(size(Model_AIC_Wavelet),find(Model_AIC_Wavelet==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_WAVE(n,Iteri).AIC=lowest_AIC;
                        GRID_WAVE(n,Iteri).FREQ=FREQ_VECTOR(B);
                        GRID_WAVE(n,Iteri).SCALING=SCALING_VECTOR(C);
                        GRID_WAVE(n,Iteri).TAPER=TAPER_VECTOR(D);
                        GRID_WAVE(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR_Wave(E);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'GaussianTaper_ModelFits_',method],'GRID_WAVE');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(FREQ_VECTOR),median(SCALING_VECTOR),median(TAPER_VECTOR),median(TIMESHIFT_VECTOR_Wave)]';
                        lower_limits = [min(FREQ_VECTOR),min(SCALING_VECTOR),min(TAPER_VECTOR),min(TIMESHIFT_VECTOR_Wave)]';
                        upper_limits = [max(FREQ_VECTOR),max(SCALING_VECTOR),max(TAPER_VECTOR),max(TIMESHIFT_VECTOR_Wave)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_WAVE(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest SSE (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_WAVE(n,Iteri).AIC=AIC_out;
                        GRAD_WAVE(n,Iteri).FREQ=params_out(1);
                        GRAD_WAVE(n,Iteri).SCALING=params_out(2);
                        GRAD_WAVE(n,Iteri).TAPER=params_out(3);
                        GRAD_WAVE(n,Iteri).TIMESHIFT=params_out(4);
                        
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'GaussianTaper_ModelFits_',method],'GRAD_WAVE');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                elseif strcmp(model(Modeli),'Gamma') && ~exist([homedir,Study,DataTitle,'Gamma_ModelFits_',method,'.mat'])
                    
                    if all(method == 'Grid')
                        tic;
                        % Gamma model
                        for SHAPEi=1:length(SHAPE_VECTOR)
                            for RATEi=1:length(RATE_VECTOR)
                                for SCALINGi=1:length(SCALING_VECTOR)
                                    for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                                        
                                        params=[SHAPE_VECTOR(SHAPEi),RATE_VECTOR(RATEi),SCALING_VECTOR(SCALINGi),TAPER_VECTOR(TAPERSHIFTi)];
                                        %                                 Model_AIC_Gamma(SHAPEi,RATEi,SCALINGi,TAPERSHIFTi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                        Model_AIC_Gamma(SHAPEi,RATEi,SCALINGi,TAPERSHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                        
                                        clear params;
                                    end % Timeshift end
                                end % Scaling end
                            end % Rate end
                        end % Shape end
                        Gammatoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Gamma);
                        n_data = dims(1).*dims(2).*dims(3).*dims(4);
                        ReshapeDaMatrix = reshape(Model_AIC_Gamma,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D,E]=ind2sub(size(Model_AIC_Gamma),find(Model_AIC_Gamma==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_GAMA(n,Iteri).AIC=lowest_AIC;
                        GRID_GAMA(n,Iteri).SHAPE=SHAPE_VECTOR(B);
                        GRID_GAMA(n,Iteri).RATE=RATE_VECTOR(C);
                        GRID_GAMA(n,Iteri).SCALE=SCALING_VECTOR(D);
                        GRID_GAMA(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR(E);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'Gamma_ModelFits_',method],'GRID_GAMA');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(SHAPE_VECTOR),median(RATE_VECTOR),median(SCALING_VECTOR),median(TIMESHIFT_VECTOR)]';
                        lower_limits = [min(SHAPE_VECTOR),min(RATE_VECTOR),min(SCALING_VECTOR),min(TIMESHIFT_VECTOR)]';
                        upper_limits = [max(SHAPE_VECTOR),max(RATE_VECTOR),max(SCALING_VECTOR),max(TIMESHIFT_VECTOR)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_GAMA(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest SSE (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_GAMA(n,Iteri).AIC=AIC_out;
                        GRAD_GAMA(n,Iteri).SHAPE=params_out(1);
                        GRAD_GAMA(n,Iteri).RATE=params_out(2);
                        GRAD_GAMA(n,Iteri).SCALE=params_out(3);
                        GRAD_GAMA(n,Iteri).TIMESHIFT=params_out(4);
                        
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'Gamma_ModelFits_',method],'GRAD_GAMA');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                elseif strcmp(model(Modeli),'TwoHalf') && ~exist([homedir,Study,DataTitle,'TwoHalfCosine_ModelFits_',method,'.mat'])
                    
                    if all(method == 'Grid')
                        tic
                        % Half Cosine model
                        for FREQi=1:length(FREQ_VECTOR)
                            for SCALINGi=1:length(SCALING_VECTOR)
                                for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                                    
                                    params=[FREQ_VECTOR(FREQi),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR(TIMESHIFTi)];
                                    % Model_AIC_Half(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},mean(Use.Pun_ERP,1),mean(Use.Rew_ERP,1),tx2disp,srate,n,Use.time1,Use.time2);
                                    Model_AIC_Half(FREQi,SCALINGi,TIMESHIFTi) = TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                                    
                                    clear params;
                                end % Timeshift end
                            end % Scaling end
                        end % Frequency end
                        Halftoc=toc;
                        
                        % -------------------------------------
                        % find "best" wave
                        dims=size(Model_AIC_Half);
                        n_data = dims(1).*dims(2).*dims(3);
                        ReshapeDaMatrix = reshape(Model_AIC_Half,1,n_data);
                        [lowest_AIC,lowest_idx]=min(ReshapeDaMatrix);
                        [B,C,D]=ind2sub(size(Model_AIC_Half),find(Model_AIC_Half==lowest_AIC,1));
                        clear dims n_data ReshapeDaMatrix
                        
                        % ########
                        GRID_2HALF(n,Iteri).AIC=lowest_AIC;
                        GRID_2HALF(n,Iteri).FREQ=FREQ_VECTOR(B);
                        GRID_2HALF(n,Iteri).SCALING=SCALING_VECTOR(C);
                        GRID_2HALF(n,Iteri).TIMESHIFT=TIMESHIFT_VECTOR(D);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'TwoHalfCosine_ModelFits_',method],'GRID_2HALF');
                        end
                        
                    elseif all(method == 'Grad')
                        tic
                        % Gradient Descent
                        
                        % for fmincon, set middle, lower, and upper values of model
                        init_params = [median(FREQ_VECTOR),median(SCALING_VECTOR),median(TIMESHIFT_VECTOR),median(FREQ_VECTOR),median(TwoModel_SCALING_VECTOR),median(TIMESHIFT_VECTOR)]';
                        lower_limits = [min(FREQ_VECTOR),min(SCALING_VECTOR),min(TIMESHIFT_VECTOR),min(FREQ_VECTOR),min(TwoModel_SCALING_VECTOR),min(TIMESHIFT_VECTOR)]';
                        upper_limits = [max(FREQ_VECTOR),max(SCALING_VECTOR),max(TIMESHIFT_VECTOR),max(FREQ_VECTOR),max(TwoModel_SCALING_VECTOR),max(TIMESHIFT_VECTOR)]';
                        
                        %
                        [xfinal,ffinal,exitflag,xstart] = rmsearch(@(params) TestRewPModel(params,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline),...
                            'fmincon',init_params,lower_limits,upper_limits,'Options',options,'InitialSample',RMRUNS);
                        GRAD_2HALF(n,Iteri).toc=toc;
                        
                        % ------------------------------------- find the values that correspond to the lowest AIC (out of RMRUNS number of starting-point iterations)
                        if sum(~isnan(ffinal))~=0
                            startidx=find(ffinal==min(ffinal));
                            params_out=xfinal(min(startidx),:);
                            AIC_out=ffinal(min(startidx),:);
                        end
                        
                        % ########
                        GRAD_2HALF(n,Iteri).AIC=AIC_out;
                        GRAD_2HALF(n,Iteri).FREQ=params_out(1);
                        GRAD_2HALF(n,Iteri).SCALING=params_out(2);
                        GRAD_2HALF(n,Iteri).TIMESHIFT=params_out(3);
                        GRAD_2HALF(n,Iteri).FREQ2=params_out(4);
                        GRAD_2HALF(n,Iteri).SCALING2=params_out(5);
                        GRAD_2HALF(n,Iteri).TIMESHIFT2=params_out(6);
                        % ########
                        
                        if n == N && Iteri == Iterations
                            save([homedir,Study,DataTitle,'TwoHalfCosine_ModelFits_',method],'GRAD_2HALF');
                        end
                        
                    else
                        fprintf('\n\nError in ''method'' variable. Redefine ''method''.\n');
                        END_SCRIPT
                    end
                    
                end % End Model testing
                
                
                
                
            end % Participant end
            
            
        end % End Model Loop
        
        clear Use
    end % End Iteration Loop
    
    clearvars -except homedir Iterations tx2disp model method RMRUNS options FREQ_VECTOR SCALING_VECTOR TIMESHIFT_VECTOR TAPER_VECTOR SHAPE_VECTOR RATE_VECTOR UseDataset
    
end % End Data Loop

%% Create the best-fitting RewP (Figures 2-5)

UseDataset = 6;

if UseDataset == 1
    % R01 Door
    Study = 'E1_SightTrain_R01Door\'; cd([homedir,Study]);
    %load('R01DoorVolt_30-Sep-2021.mat');
    load('R01DepDoorVolt_15-Aug-2022.mat');
    srate = 500;
    DataTitle = 'R01_Door_';
    CondiTitle = {'Sight'};
    Iterations = 1;
    LossCondi = 3; % Specific to this Doors Task. 2 = Neutral; 3 = Loss;
    Electrode = 36;
    N = size(DoorVolt.ERPs_FB,1);
    
    % Pull info for model
    
    % Individual differences
    Training_Sight.Pun_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,LossCondi) );
    Training_Sight.Rew_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,1) );
    
    % Everything Else
    Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
    Training_Sight.time1 = 100;
    Training_Sight.time2 = 500;
    % Training_Sight.BDI = DoorVolt.BDI';
    % Training_Sight.BIS = DoorVolt.BIS';
    % Training_Sight.BAS_D = DoorVolt.BAS_D';
    % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sight.age = DoorVolt.Age';
    % Training_Sight.sex = DoorVolt.sex';
    
elseif UseDataset == 2
    % Sight Or Sound
    Study = 'E2_SightOrSound\'; cd([homedir,Study]);
    load('SoS_Data.mat','MEGA_ERPs');
    srate = 500;
    DataTitle = 'SightOrSound_';
    CondiTitle = {'Sight','Sound'};
    Iterations = 2;
    LossCondi = [4 1];
    N = size(MEGA_ERPs.FBbyCUE,1);
    Electrode = 36;
    
    % Pull info for model
    
    % Individual differences
    Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi(1)) );
    Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2) );
    
    % Everything Else
    Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
    Training_Sight.time1 = 100;
    Training_Sight.time2 = 500;
    Training_Sight.UseBaseline = 0;
    % Training_Sight.BDI = DoorVolt.BDI';
    % Training_Sight.BIS = DoorVolt.BIS';
    % Training_Sight.BAS_D = DoorVolt.BAS_D';
    % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sight.age = DoorVolt.Age';
    % Training_Sight.sex = DoorVolt.sex';
    
    % Individual differences
    Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi(2)) );
    Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3) );
    
    % Everything Else
    Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
    Training_Sound.time1 = 100;
    Training_Sound.time2 = 500;
    % Training_Sound.BDI = DoorVolt.BDI';
    % Training_Sound.BIS = DoorVolt.BIS';
    % Training_Sound.BAS_D = DoorVolt.BAS_D';
    % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sound.age = DoorVolt.Age';
    % Training_Sound.sex = DoorVolt.sex';
    
elseif UseDataset == 3
    
    % Sight And Sound
    Study = 'E3_SightAndSound\'; cd([homedir,Study]);
    load('SnS_Data.mat');
    srate = 500;
    DataTitle = 'SightAndSound_';
    CondiTitle = {'Sight','Sound','Combo'};
    Iterations = 3;
    LossCondi = 1;
    N = size(MEGA_ERPs.FBbyCUE,1);
    Electrode = 36;
    
    % Pull info for model
    
    % Individual differences
    Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
    Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3) );
    
    % Everything Else
    Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
    Training_Sight.time1 = 100;
    Training_Sight.time2 = 500;
    Training_Sight.UseBaseline = 0;
    % Training_Sight.BDI = DoorVolt.BDI';
    % Training_Sight.BIS = DoorVolt.BIS';
    % Training_Sight.BAS_D = DoorVolt.BAS_D';
    % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sight.age = DoorVolt.Age';
    % Training_Sight.sex = DoorVolt.sex';
    
    % Individual differences
    Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
    Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,4) );
    
    % Everything Else
    Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
    Training_Sound.time1 = 100;
    Training_Sound.time2 = 500;
    Training_Sound.UseBaseline = 0;
    % Training_Sound.BDI = DoorVolt.BDI';
    % Training_Sound.BIS = DoorVolt.BIS';
    % Training_Sound.BAS_D = DoorVolt.BAS_D';
    % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sound.age = DoorVolt.Age';
    % Training_Sound.sex = DoorVolt.sex';
    
    % Individual differences
    Training_Combo.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi) );
    Training_Combo.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2) );
    
    % Everything Else
    Training_Combo.RewP_Diff = Training_Combo.Rew_ERP - Training_Combo.Pun_ERP;
    Training_Combo.time1 = 100;
    Training_Combo.time2 = 500;
    Training_Combo.UseBaseline = 0;
    % Training_Combo.BDI = DoorVolt.BDI';
    % Training_Combo.BIS = DoorVolt.BIS';
    % Training_Combo.BAS_D = DoorVolt.BAS_D';
    % Training_Combo.BAS_FS = DoorVolt.BAS_FS';
    % Training_Combo.BAS_RR = DoorVolt.BAS_RR';
    % Training_Combo.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Combo.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Combo.age = DoorVolt.Age';
    % Training_Combo.sex = DoorVolt.sex';
    
elseif UseDataset==4
    
    % Gremlin
    Study = 'E4_SoundDelay_Gremlin\'; cd([homedir,Study]);
    load('V_Gremlin_16-Aug-2022.mat');
    srate = 500;
    DataTitle = 'V_Gremlin_16-Aug-2022';
    CondiTitle = {'100 ms','300 ms','500 ms'};
    Iterations = 3;
    LossCondi = 1;
    N = size(MEGA_ERPs.FBbyCUE,1);
    Electrode = 21;
    
    % Pull info for model
    
    % Individual differences
    Training_Sight.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,1) );
    Training_Sight.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,1) );
    
    % Everything Else
    Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
    Training_Sight.time1 = 100;
    Training_Sight.time2 = 600;
    Training_Sight.UseBaseline = 100;
    % Training_Sight.BDI = DoorVolt.BDI';
    % Training_Sight.BIS = DoorVolt.BIS';
    % Training_Sight.BAS_D = DoorVolt.BAS_D';
    % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sight.age = DoorVolt.Age';
    % Training_Sight.sex = DoorVolt.sex';
    
    % Individual differences
    Training_Sound.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,2) );
    Training_Sound.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,2) );
    
    % Everything Else
    Training_Sound.RewP_Diff = Training_Sound.Rew_ERP - Training_Sound.Pun_ERP;
    Training_Sound.time1 = 300;
    Training_Sound.time2 = 800;
    Training_Sound.UseBaseline = 300;
    % Training_Sound.BDI = DoorVolt.BDI';
    % Training_Sound.BIS = DoorVolt.BIS';
    % Training_Sound.BAS_D = DoorVolt.BAS_D';
    % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sound.age = DoorVolt.Age';
    % Training_Sound.sex = DoorVolt.sex';
    
    % Individual differences
    Training_Combo.Pun_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,LossCondi,3) );
    Training_Combo.Rew_ERP = squeeze( MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,3) );
    
    % Everything Else
    Training_Combo.RewP_Diff = Training_Combo.Rew_ERP - Training_Combo.Pun_ERP;
    Training_Combo.time1 = 500;
    Training_Combo.time2 = 1000;
    Training_Combo.UseBaseline = 500;
    % Training_Sound.BDI = DoorVolt.BDI';
    % Training_Sound.BIS = DoorVolt.BIS';
    % Training_Sound.BAS_D = DoorVolt.BAS_D';
    % Training_Sound.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sound.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sound.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sound.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sound.age = DoorVolt.Age';
    % Training_Sound.sex = DoorVolt.sex';
    
elseif UseDataset == 5
    % R01 Depressed Door
    Study = 'T1_SightDep_R01Door\'; cd([homedir,Study]);
    %load('R01DoorVolt_30-Sep-2021.mat');
    load('R01DepDoorVolt_14-Mar-2023.mat');
    srate = 500;
    DataTitle = 'R01_Door_';
    CondiTitle = {'Sight'};
    Iterations = 1;
    LossCondi = 3; % Specific to this Doors Task. 2 = Neutral; 3 = Loss;
    Electrode = 36;
    N = size(DoorVolt.ERPs_FB,1);
    
    % Pull info for model
    
    % Individual differences
    Training_Sight.Pun_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,LossCondi) );
    Training_Sight.Rew_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,1:751,1) );
    
    % Everything Else
    Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
    Training_Sight.time1 = 100;
    Training_Sight.time2 = 500;
    Training_Sight.UseBaseline = 0;
    % Training_Sight.BDI = DoorVolt.BDI';
    % Training_Sight.BIS = DoorVolt.BIS';
    % Training_Sight.BAS_D = DoorVolt.BAS_D';
    % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
    % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
    % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
    % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
    % Training_Sight.age = DoorVolt.Age';
    % Training_Sight.sex = DoorVolt.sex';
    
    elseif UseDataset == 6
        
        % Cavanagh et al., 2019
        Study = 'T2_RL_Dep\'; cd([homedir,Study]);
        %load('R01DoorVolt_30-Sep-2021.mat');
        load('PS_RewSens.mat');
        srate = 500;
        DataTitle = 'PST_Depression_Cav2019_';
        CondiTitle = {'Sight'};
        Iterations = 1;
        LossCondi = 1; 
        WinCondi=2;
        Electrode = 36;
        DepressedParticipants=DEP;
        N = size(MEGA_ERP,1);
        
        % Pull info for model
        
        % Individual differences
        Training_Sight.Pun_ERP = squeeze( MEGA_ERP(:,Electrode,1:751,LossCondi) );
        Training_Sight.Rew_ERP = squeeze( MEGA_ERP(:,Electrode,1:751,WinCondi) );
        
        % Everything Else
        Training_Sight.RewP_Diff = Training_Sight.Rew_ERP - Training_Sight.Pun_ERP;
        Training_Sight.time1 = 100;
        Training_Sight.time2 = 500;
        Training_Sight.UseBaseline = 0;
        % Training_Sight.BDI = DoorVolt.BDI';
        % Training_Sight.BIS = DoorVolt.BIS';
        % Training_Sight.BAS_D = DoorVolt.BAS_D';
        % Training_Sight.BAS_FS = DoorVolt.BAS_FS';
        % Training_Sight.BAS_RR = DoorVolt.BAS_RR';
        % Training_Sight.MASQ_AD = DoorVolt.MASQ_AD';
        % Training_Sight.MASQ_GDD = DoorVolt.MASQ_GDD';
        % Training_Sight.age = DoorVolt.Age';
        % Training_Sight.sex = DoorVolt.sex';
    
end

PlotAIC = zeros(size(Training_Sight.Pun_ERP,1),size(model,2),Iterations);

for Modaliti = 1:Iterations
    
    % Pick one to use
    if Modaliti == 1
        Use = Training_Sight;
    elseif Modaliti == 2
        Use = Training_Sound;
    elseif Modaliti == 3
        Use = Training_Combo;
    end 
    
    Use.AvgPun = mean(Use.Pun_ERP,1);
    Use.AvgRew = mean(Use.Rew_ERP,1);
    Use.AvgDif = mean(Use.RewP_Diff,1);
    ylimit = [-12 20];
    
    figure; hold on
    plot(tx2disp,Use.AvgPun,'r','linewidth',2);
    plot(tx2disp,Use.AvgRew,'g','linewidth',2);
    plot(tx2disp,Use.AvgDif,'b','linewidth',2);
    set(gca,'ylim',ylimit)
    title(['Real Waveforms']);
    
    figure;
    
    % Create Artificial bursts and constructed waveforms for plotting
    for Modeli = 1:size(model,2)
        
        if Modeli == 1
            PlotSpace = [1 2 3 4];
        elseif Modeli > 1
            PlotSpace = [1+4.*(Modeli-1) 2+4.*(Modeli-1) 3+4.*(Modeli-1) 4+4.*(Modeli-1)];
        end
        
        % Comment in to use average
%         n=1;
        
        if strcmp(model(Modeli),'Full')
            
            load([homedir,Study,DataTitle,'FullCosine_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_FULL = GRAD_FULL;
                clear GRAD_FULL
            elseif all(method=='Grid')
                MODEL_FULL = GRID_FULL;
                clear GRID_FULL
            end
            
            NumOfParams = 3;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.Full(parti,1)=MODEL_FULL(parti,Modaliti).AIC;
                SPSS{Modaliti}.Full(parti,2)=MODEL_FULL(parti,Modaliti).FREQ;
                SPSS{Modaliti}.Full(parti,3)=MODEL_FULL(parti,Modaliti).SCALING;
                SPSS{Modaliti}.Full(parti,4)=MODEL_FULL(parti,Modaliti).TIMESHIFT;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_FULL(n,Modaliti).FREQ MODEL_FULL(n,Modaliti).SCALING MODEL_FULL(n,Modaliti).TIMESHIFT];
                [FULL_ArtificialBurst(Modaliti,n,:), FULL_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_FULL(n,Modaliti).AIC]';
                SPSS{Modaliti}.Full(n,5)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            FULL_params(Modaliti,:) = [mean([MODEL_FULL(:,Modaliti).FREQ],2) mean([MODEL_FULL(:,Modaliti).SCALING],2) round(mean([MODEL_FULL(:,Modaliti).TIMESHIFT],2))];
%             [FULL_ArtificialBurst(Modaliti,n,:), FULL_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(FULL_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
%             
%             PlotAIC(:,Modeli,Modaliti) = [MODEL_FULL(:,Modaliti).AIC]';
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(FULL_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_FULL(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['FREQ: ',num2str( FULL_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['SCAL: ',num2str( FULL_params(Modaliti,2) )],'sc');
            text(.1,.35, ['TIME: ',num2str( FULL_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(FULL_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(FULL_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' FULL Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(FULL_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' FULL Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(FULL_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(FULL_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear FULL_ArtificialBurst FULL_ConstructedWaveform FULL_params MODEL_FULL ModelAIC
            
        elseif strcmp(model(Modeli),'Half')
            
            load([homedir,Study,DataTitle,'HalfCosine_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_HALF = GRAD_HALF;
                clear GRAD_HALF
            elseif all(method=='Grid')
                MODEL_HALF = GRID_HALF;
                clear GRID_HALF
            end
                  
            NumOfParams = 3;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.Half(parti,1)=MODEL_HALF(parti,Modaliti).AIC;
                SPSS{Modaliti}.Half(parti,2)=MODEL_HALF(parti,Modaliti).FREQ;
                SPSS{Modaliti}.Half(parti,3)=MODEL_HALF(parti,Modaliti).SCALING;
                SPSS{Modaliti}.Half(parti,4)=MODEL_HALF(parti,Modaliti).TIMESHIFT;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_HALF(n,Modaliti).FREQ MODEL_HALF(n,Modaliti).SCALING MODEL_HALF(n,Modaliti).TIMESHIFT];
                [HALF_ArtificialBurst(Modaliti,n,:), HALF_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_HALF(n,Modaliti).AIC]';
                SPSS{Modaliti}.Half(n,5)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            HALF_params(Modaliti,:) = [mean([MODEL_HALF(:,Modaliti).FREQ],2) mean([MODEL_HALF(:,Modaliti).SCALING],2) round(mean([MODEL_HALF(:,Modaliti).TIMESHIFT],2))];
%             [HALF_ArtificialBurst(Modaliti,n,:), HALF_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(HALF_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(HALF_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_HALF(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['FREQ: ',num2str( HALF_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['SCAL: ',num2str( HALF_params(Modaliti,2) )],'sc');
            text(.1,.35, ['TIME: ',num2str( HALF_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' HALF Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' HALF Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(HALF_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(HALF_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear HALF_ArtificialBurst HALF_ConstructedWaveform HALF_params MODEL_HALF ModelAIC
            
        elseif strcmp(model(Modeli),'Gauss')
            
            load([homedir,Study,DataTitle,'Gaussian_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_GAUS = GRAD_GAUS;
                clear GRAD_GAUS
            elseif all(method=='Grid')
                MODEL_GAUS = GRID_GAUS;
                clear GRID_GAUS
            end
            
            NumOfParams = 3;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.Gauss(parti,1)=MODEL_GAUS(parti,Modaliti).AIC;
                SPSS{Modaliti}.Gauss(parti,2)=MODEL_GAUS(parti,Modaliti).FREQ;
                SPSS{Modaliti}.Gauss(parti,3)=MODEL_GAUS(parti,Modaliti).SCALING;
                SPSS{Modaliti}.Gauss(parti,4)=MODEL_GAUS(parti,Modaliti).TIMESHIFT;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_GAUS(n,Modaliti).FREQ MODEL_GAUS(n,Modaliti).SCALING MODEL_GAUS(n,Modaliti).TIMESHIFT];
                [GAUS_ArtificialBurst(Modaliti,n,:), GAUS_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_GAUS(n,Modaliti).AIC]';
                SPSS{Modaliti}.Gauss(n,5)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            GAUS_params(Modaliti,:) = [mean([MODEL_GAUS(:,Modaliti).FREQ],2) mean([MODEL_GAUS(:,Modaliti).SCALING],2) round(mean([MODEL_GAUS(:,Modaliti).TIMESHIFT],2))];
%             [GAUS_ArtificialBurst(Modaliti,n,:), GAUS_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(GAUS_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(GAUS_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_GAUS(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['FREQ: ',num2str( GAUS_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['SCAL: ',num2str( GAUS_params(Modaliti,2) )],'sc');
            text(.1,.35, ['TIME: ',num2str( GAUS_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(GAUS_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(GAUS_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' GAUSS Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(GAUS_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' GAUSS Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(GAUS_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(GAUS_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear GAUS_ArtificialBurst GAUS_ConstructedWaveform GAUS_params MODEL_GAUS ModelAIC
            
        elseif strcmp(model(Modeli),'Wavelet')
            
            load([homedir,Study,DataTitle,'GaussianTaper_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_WAVE = GRAD_WAVE;
                clear GRAD_WAVE
            elseif all(method=='Grid')
                MODEL_WAVE = GRID_WAVE;
                clear GRID_WAVE
            end
            
            NumOfParams = 4;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.Wavelet(parti,1)=MODEL_WAVE(parti,Modaliti).AIC;
                SPSS{Modaliti}.Wavelet(parti,2)=MODEL_WAVE(parti,Modaliti).FREQ;
                SPSS{Modaliti}.Wavelet(parti,3)=MODEL_WAVE(parti,Modaliti).SCALING;
                SPSS{Modaliti}.Wavelet(parti,4)=MODEL_WAVE(parti,Modaliti).TIMESHIFT;
                SPSS{Modaliti}.Wavelet(parti,5)=MODEL_WAVE(parti,Modaliti).TAPER;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_WAVE(n,Modaliti).FREQ MODEL_WAVE(n,Modaliti).SCALING MODEL_WAVE(n,Modaliti).TAPER MODEL_WAVE(n,Modaliti).TIMESHIFT];
                [WAVE_ArtificialBurst(Modaliti,n,:), WAVE_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_WAVE(n,Modaliti).AIC]';
                SPSS{Modaliti}.Wavelet(n,6)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            WAVE_params(Modaliti,:) = [mean([MODEL_WAVE(:,Modaliti).FREQ],2) mean([MODEL_WAVE(:,Modaliti).SCALING],2) round(mean([MODEL_WAVE(:,Modaliti).TAPER],2)) round(mean([MODEL_WAVE(:,Modaliti).TIMESHIFT],2))];
%             [WAVE_ArtificialBurst(Modaliti,n,:), WAVE_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(WAVE_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(WAVE_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_WAVE(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['FREQ: ',num2str( WAVE_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['SCAL: ',num2str( WAVE_params(Modaliti,2) )],'sc');
            text(.1,.35, ['TAPR: ',num2str( WAVE_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['TIME: ',num2str( WAVE_params(Modaliti,4)-100 )],'sc');
            text(.1,.05,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(WAVE_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(WAVE_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' WAVELET Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(WAVE_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' WAVELET Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(WAVE_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(WAVE_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear WAVE_ArtificialBurst WAVE_ConstructedWaveform WAVE_params MODEL_WAVE ModelAIC
            
        elseif strcmp(model(Modeli),'Gamma')
            
            load([homedir,Study,DataTitle,'Gamma_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_GAMA = GRAD_GAMA;
                clear GRAD_GAMA
            elseif all(method=='Grid')
                MODEL_GAMA = GRID_GAMA;
                clear GRID_GAMA
            end
            
            NumOfParams = 4;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.Gamma(parti,1)=MODEL_GAMA(parti,Modaliti).AIC;
                SPSS{Modaliti}.Gamma(parti,2)=MODEL_GAMA(parti,Modaliti).SHAPE;
                SPSS{Modaliti}.Gamma(parti,3)=MODEL_GAMA(parti,Modaliti).RATE;
                SPSS{Modaliti}.Gamma(parti,4)=MODEL_GAMA(parti,Modaliti).SCALE;
                SPSS{Modaliti}.Gamma(parti,5)=MODEL_GAMA(parti,Modaliti).TIMESHIFT;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_GAMA(n,Modaliti).SHAPE MODEL_GAMA(n,Modaliti).RATE MODEL_GAMA(n,Modaliti).SCALE MODEL_GAMA(n,Modaliti).TIMESHIFT];
                [GAMA_ArtificialBurst(Modaliti,n,:) GAMA_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_GAMA(n,Modaliti).AIC]';
                SPSS{Modaliti}.Gamma(n,6)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            GAMA_params(Modaliti,:) = [mean([MODEL_GAMA(:,Modaliti).SHAPE],2) mean([MODEL_GAMA(:,Modaliti).RATE],2) mean([MODEL_GAMA(:,Modaliti).SCALE],2) round(mean([MODEL_GAMA(:,Modaliti).TIMESHIFT],2))];
%             [GAMA_ArtificialBurst(Modaliti,n,:) GAMA_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(GAMA_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(GAMA_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_GAMA(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['SHAP: ',num2str( GAMA_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['RATE: ',num2str( GAMA_params(Modaliti,2) )],'sc');
            text(.1,.35, ['SCAL: ',num2str( GAMA_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['TIME: ',num2str( GAMA_params(Modaliti,4) )],'sc');
            text(.1,.05,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(GAMA_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(GAMA_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' GAMMA Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(GAMA_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' GAMMA Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(GAMA_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(GAMA_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear GAMA_ArtificialBurst GAMA_ConstructedWaveform GAMA_params MODEL_GAMA ModelAIC
            
        elseif strcmp(model(Modeli),'TwoHalf')
            
            load([homedir,Study,DataTitle,'TwoHalfCosine_ModelFits_',method])
            
            if all(method=='Grad')
                MODEL_2HALF = GRAD_2HALF;
                clear GRAD_2HALF
            elseif all(method=='Grid')
                MODEL_2HALF = GRID_2HALF;
                clear GRID_2HALF
            end
                  
            NumOfParams = 3;
            NumOfSamples = (Use.time2-Use.time1)./2;
            for parti = 1:N
                SPSS{Modaliti}.TwoHalf(parti,1)=MODEL_2HALF(parti,Modaliti).AIC;
                SPSS{Modaliti}.TwoHalf(parti,2)=MODEL_2HALF(parti,Modaliti).FREQ;
                SPSS{Modaliti}.TwoHalf(parti,3)=MODEL_2HALF(parti,Modaliti).SCALING;
                SPSS{Modaliti}.TwoHalf(parti,4)=MODEL_2HALF(parti,Modaliti).TIMESHIFT;
                SPSS{Modaliti}.TwoHalf(parti,5)=MODEL_2HALF(parti,Modaliti).FREQ2;
                SPSS{Modaliti}.TwoHalf(parti,6)=MODEL_2HALF(parti,Modaliti).SCALING2;
                SPSS{Modaliti}.TwoHalf(parti,7)=MODEL_2HALF(parti,Modaliti).TIMESHIFT2;
            end
            
            % Average the added waveform
            for n=1:N
                ParticipantParams = [MODEL_2HALF(n,Modaliti).FREQ MODEL_2HALF(n,Modaliti).SCALING MODEL_2HALF(n,Modaliti).TIMESHIFT MODEL_2HALF(n,Modaliti).FREQ2 MODEL_2HALF(n,Modaliti).SCALING2 MODEL_2HALF(n,Modaliti).TIMESHIFT2];
                [HALF2_ArtificialBurst(Modaliti,n,:), HALF2_ConstructedWaveform(Modaliti,n,:), AIC(n,:), BurstMaximum] = Recreate_TestRewPModel(ParticipantParams,model{Modeli},Use.Pun_ERP,Use.Rew_ERP,tx2disp,srate,n,Use.time1,Use.time2,Use.UseBaseline);
                
                PlotAIC(n,Modeli,Modaliti) = [MODEL_2HALF(n,Modaliti).AIC]';
                SPSS{Modaliti}.Half(n,5)=BurstMaximum;
                
                clear ParticipantParams BurstMaximum
            end
            
            % Average participant params and create waveform
            HALF2_params(Modaliti,:) = [mean([MODEL_2HALF(:,Modaliti).FREQ],2) mean([MODEL_2HALF(:,Modaliti).SCALING],2) round(mean([MODEL_2HALF(:,Modaliti).TIMESHIFT],2)) mean([MODEL_2HALF(:,Modaliti).FREQ2],2) mean([MODEL_2HALF(:,Modaliti).SCALING2],2) round(mean([MODEL_2HALF(:,Modaliti).TIMESHIFT2],2))];
%             [HALF_ArtificialBurst(Modaliti,n,:), HALF_ConstructedWaveform(Modaliti,n,:) ModelAIC] = Recreate_TestRewPModel(HALF_params(Modaliti,:),model{Modeli},Use.AvgPun,Use.AvgRew,tx2disp,srate,n,Use.time1,Use.time2);
            
            T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
            SSE_Full = sum((Use.AvgRew(:,T1:T2)-squeeze(mean(HALF2_ConstructedWaveform(Modaliti,:,T1:T2),2))').^2);
            ModelAIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;
            clear NumOfSamples NumOfParams SSE_Full T1 T2 
            
            % Params
            subplot(size(model,2),4,PlotSpace(1)); hold on;
            text(.1,.8,  ['AIC: ', num2str(mean([MODEL_2HALF(:,Modaliti).AIC],2))],'sc');
            text(.1,.65, ['FREQ: ',num2str( HALF2_params(Modaliti,1) )],'sc');
            text(.1,.5,  ['SCAL: ',num2str( HALF2_params(Modaliti,2) )],'sc');
            text(.1,.35, ['TIME: ',num2str( HALF2_params(Modaliti,3) )],'sc');
            text(.1,.2,  ['Actual AIC: ',num2str( ModelAIC )],'sc');
            
            % Artificial Waveforms
            subplot(size(model,2),4,PlotSpace(2)); hold on;
            plot(tx2disp,Use.AvgPun,'r','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF2_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF2_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' HALF Artificial Waveforms and Punishment']);
            
            % RewP Comparison
            subplot(size(model,2),4,PlotSpace(3)); hold on;
            rectangle('Position',[Use.time1,-10,(Use.time2-Use.time1),30],'Curvature',0.1,'FaceColor',[.5 .5 .5])
            plot(tx2disp,Use.AvgRew,'g','linewidth',2);
            plot(tx2disp,squeeze(mean(HALF2_ConstructedWaveform(Modaliti,:,:),2)),'k','linewidth',2);
            set(gca,'ylim',ylimit)
            title([CondiTitle{Modaliti},' HALF Artificial Waveforms and Reward']);
            
            % Subject Comparison
            subplot(size(model,2),4,PlotSpace(4)); hold on;
            for  n=1:N
                plot(tx2disp,squeeze(HALF2_ArtificialBurst(Modaliti,n,:)),':','linewidth',1.5);
            end
            plot(tx2disp,squeeze(mean(HALF2_ArtificialBurst(Modaliti,:,:),2)),'c','linewidth',2);
            
            clear HALF2_ArtificialBurst HALF2_ConstructedWaveform HALF2_params MODEL_2HALF ModelAIC
            
        end
        
    end
end

%% Bayesian Model Selection

for Modaliti = 1:Iterations

    X1 = squeeze(mean(PlotAIC(:,1,Modaliti),1)); X2 = squeeze(mean(PlotAIC(:,2,Modaliti),1)); X3 = squeeze(mean(PlotAIC(:,3,Modaliti),1)); ...
        X4 = squeeze(mean(PlotAIC(:,4,Modaliti),1)); X5 = squeeze(mean(PlotAIC(:,5,Modaliti),1));
    
    AverageAIC = [X1,X2,X3,X4,X5];
    
    % Find best fitting model
    ModelIndexing = 1:5;
    [BestModelAIC, BestModelLoc] = min(AverageAIC);
    if Modaliti == 2
        BestModelLoc = 1;
    end
    IndexOtherModels = ModelIndexing(ModelIndexing~=BestModelLoc);
    
    figure;
    hold on;
    plot(AverageAIC,'bs');
    title(['AIC by Model, condition = ',num2str(Modaliti)]);
    
    % Scatter for AIC
    for ai=1:5
        plot(ai-.15:.30/(length(PlotAIC)-1):ai+.15,PlotAIC(:,ai),'m.');
    end
    
    % Bayesian Model Selection
    for xxi=1:size(IndexOtherModels,2);
        TestModel = IndexOtherModels(xxi);
        [alpha, exp_r, xp] = spm_BMS(squeeze(PlotAIC(:,[BestModelLoc,TestModel],Modaliti)), [], 1, 1, 1);
        BMS(xxi,:)=[alpha, exp_r, xp];
        clear alpha exp_r xp TestModel;
    end
    
    for i = 1:size(PlotAIC,1)
        RankAIC(i,:) = floor(tiedrank([squeeze(PlotAIC(i,:,Modaliti))]));
    end
    
    FullWin = sum(RankAIC(:,1)==1)
    HalfWin = sum(RankAIC(:,2)==1)
    GausWin = sum(RankAIC(:,3)==1)
    WaveWin = sum(RankAIC(:,4)==1)
    GamaWin = sum(RankAIC(:,5)==1)
    
    clear X1 X2 X3 X4 X5 ModelIndexing BestModelAIC BestModelLoc AverageAIC IndexOtherModels BMS RankAIC FullWin HalfWin GausWin WaveWin GamaWin
end

%% Tables
MeanTable = zeros(size(model,2),7,Iterations); SDTable = zeros(size(model,2),7,Iterations);
for Modaliti = 1:Iterations
    
    % AIC
    MeanTable(1,1,Modaliti) = mean(SPSS{1,Modaliti}.Full(:,1),1);    SDTable(1,1,Modaliti) = std(SPSS{1,Modaliti}.Full(:,1),[],1);
    MeanTable(2,1,Modaliti) = mean(SPSS{1,Modaliti}.Half(:,1),1);    SDTable(2,1,Modaliti) = std(SPSS{1,Modaliti}.Half(:,1),[],1);
    MeanTable(3,1,Modaliti) = mean(SPSS{1,Modaliti}.Gauss(:,1),1);   SDTable(3,1,Modaliti) = std(SPSS{1,Modaliti}.Gauss(:,1),[],1);
    MeanTable(4,1,Modaliti) = mean(SPSS{1,Modaliti}.Wavelet(:,1),1); SDTable(4,1,Modaliti) = std(SPSS{1,Modaliti}.Wavelet(:,1),[],1);
    MeanTable(5,1,Modaliti) = mean(SPSS{1,Modaliti}.Gamma(:,1),1);   SDTable(5,1,Modaliti) = std(SPSS{1,Modaliti}.Gamma(:,1),[],1);
    
    % Frequency
    MeanTable(1,2,Modaliti) = mean(SPSS{1,Modaliti}.Full(:,2),1);    SDTable(1,2,Modaliti) = std(SPSS{1,Modaliti}.Full(:,2),[],1);
    MeanTable(2,2,Modaliti) = mean(SPSS{1,Modaliti}.Half(:,2),1);    SDTable(2,2,Modaliti) = std(SPSS{1,Modaliti}.Half(:,2),[],1);
    MeanTable(3,2,Modaliti) = mean(SPSS{1,Modaliti}.Gauss(:,2),1);   SDTable(3,2,Modaliti) = std(SPSS{1,Modaliti}.Gauss(:,2),[],1);
    MeanTable(4,2,Modaliti) = mean(SPSS{1,Modaliti}.Wavelet(:,2),1); SDTable(4,2,Modaliti) = std(SPSS{1,Modaliti}.Wavelet(:,2),[],1);
    MeanTable(5,2,Modaliti) = 0;                                     SDTable(5,2,Modaliti) = 0;
    
    % Shape
    MeanTable(1,3,Modaliti) = 0;                                     SDTable(1,3,Modaliti) = 0;
    MeanTable(2,3,Modaliti) = 0;                                     SDTable(2,3,Modaliti) = 0;
    MeanTable(3,3,Modaliti) = 0;                                     SDTable(3,3,Modaliti) = 0;
    MeanTable(4,3,Modaliti) = 0;                                     SDTable(4,3,Modaliti) = 0;
    MeanTable(5,3,Modaliti) = mean(SPSS{1,Modaliti}.Gamma(:,2),1);   SDTable(5,3,Modaliti) = std(SPSS{1,Modaliti}.Gamma(:,2),[],1);
    
    % Rate
    MeanTable(1,4,Modaliti) = 0;                                     SDTable(1,4,Modaliti) = 0;
    MeanTable(2,4,Modaliti) = 0;                                     SDTable(2,4,Modaliti) = 0;
    MeanTable(3,4,Modaliti) = 0;                                     SDTable(3,4,Modaliti) = 0;
    MeanTable(4,4,Modaliti) = 0;                                     SDTable(4,4,Modaliti) = 0;
    MeanTable(5,4,Modaliti) = mean(SPSS{1,Modaliti}.Gamma(:,3),1);   SDTable(5,4,Modaliti) = std(SPSS{1,Modaliti}.Gamma(:,3),[],1);
    
    % Taper
    MeanTable(1,5,Modaliti) = 0;                                     SDTable(1,5,Modaliti) = 0;
    MeanTable(2,5,Modaliti) = 0;                                     SDTable(2,5,Modaliti) = 0;
    MeanTable(3,5,Modaliti) = 0;                                     SDTable(3,5,Modaliti) = 0;
    MeanTable(4,5,Modaliti) = mean(SPSS{1,Modaliti}.Wavelet(:,5),1); SDTable(4,5,Modaliti) = std(SPSS{1,Modaliti}.Wavelet(:,5),[],1);
%     MeanTable(5,5,Modaliti) = 0;                                     SDTable(5,5,Modaliti) = 0;
    
    % Scaling
    MeanTable(1,6,Modaliti) = mean(SPSS{1,Modaliti}.Full(:,3),1);    SDTable(1,6,Modaliti) = std(SPSS{1,Modaliti}.Full(:,3),[],1);
    MeanTable(2,6,Modaliti) = mean(SPSS{1,Modaliti}.Half(:,3),1);    SDTable(2,6,Modaliti) = std(SPSS{1,Modaliti}.Half(:,3),[],1);
    MeanTable(3,6,Modaliti) = mean(SPSS{1,Modaliti}.Gauss(:,3),1);   SDTable(3,6,Modaliti) = std(SPSS{1,Modaliti}.Gauss(:,3),[],1);
    MeanTable(4,6,Modaliti) = mean(SPSS{1,Modaliti}.Wavelet(:,3),1); SDTable(4,6,Modaliti) = std(SPSS{1,Modaliti}.Wavelet(:,3),[],1);
    MeanTable(5,6,Modaliti) = mean(SPSS{1,Modaliti}.Gamma(:,4),1);   SDTable(5,6,Modaliti) = std(SPSS{1,Modaliti}.Gamma(:,4),[],1);
    
    % Burst Maximum
    MeanTable(1,7,Modaliti) = mean(SPSS{1,Modaliti}.Full(:,5),1);    SDTable(1,7,Modaliti) = std(SPSS{1,Modaliti}.Full(:,5),[],1);
    MeanTable(2,7,Modaliti) = mean(SPSS{1,Modaliti}.Half(:,5),1);    SDTable(2,7,Modaliti) = std(SPSS{1,Modaliti}.Half(:,5),[],1);
    MeanTable(3,7,Modaliti) = mean(SPSS{1,Modaliti}.Gauss(:,5),1);   SDTable(3,7,Modaliti) = std(SPSS{1,Modaliti}.Gauss(:,5),[],1);
    MeanTable(4,7,Modaliti) = mean(SPSS{1,Modaliti}.Wavelet(:,6),1); SDTable(4,7,Modaliti) = std(SPSS{1,Modaliti}.Wavelet(:,6),[],1);
    MeanTable(5,7,Modaliti) = mean(SPSS{1,Modaliti}.Gamma(:,6),1);   SDTable(5,7,Modaliti) = std(SPSS{1,Modaliti}.Gamma(:,6),[],1);
    
end



%% Figure 1. Sample Models

% ~~~~~~~~~~ Base Params for figure ~~~~~~~~~~ %
BaseFrex = 2;
BaseScal = 10;
BaseTime = 100;
BaseTapr = 45;
BaseShpe = 7;
BaseRate = 0.1;
srate = 500;
ylimit = [-10 20];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% SHAPE_VECTOR=           5:5:40;
% RATE_VECTOR=            0.2:0.05:1;

figure;

for Modeli = 1:size(model,2)
        
        if Modeli == 1
            PlotSpace = [1 2 3];
        elseif Modeli > 1
            PlotSpace = [1+3.*(Modeli-1) 2+3.*(Modeli-1) 3+3.*(Modeli-1)];
        end
        
        if strcmp(model(Modeli),'Full')
            
            % Get time
            full_ti= -((1/BaseFrex)*.25):(1000/srate)/1000:((1/BaseFrex)*.75);
            
            % Fit Freq
            full_burst = cos(2*pi*BaseFrex.*full_ti);
            
            % Amplitude scaling
            SCALE_Full = full_burst.*BaseScal;
            
            % Add Baseline and remove tail if beyond time window
            AddBaseline  = 251+BaseTime;
            TIME_Full = [zeros(1,AddBaseline),SCALE_Full,zeros(1,(751-(length(SCALE_Full)+AddBaseline)))];
            TIME_Full = TIME_Full(:,1:751);
            
            % Create Artificial Burst and Constructed Waveform for plotting
            ArtificialBurst = TIME_Full;
            
            ModelTitle = 'Full Cosine';
            
            clear full_ti full_burst SCALE_Full AddBaseline TIME_Full
            
        elseif strcmp(model(Modeli),'Half')
            
            % Get time
            half_ti= -((1/BaseFrex)*.25):(1000/srate)/1000:((1/BaseFrex)*.25);
            
            % Fit Freq
            half_burst = cos(2*pi*BaseFrex.*half_ti);
            
            % Amplitude scaling
            SCALE_Half = half_burst.*BaseScal;
            
            % Add Baseline and remove tail if beyond time window
            AddBaseline  = 251+BaseTime;
            TIME_Half = [zeros(1,AddBaseline),SCALE_Half,zeros(1,(751-(length(SCALE_Half)+AddBaseline)))];
            TIME_Half = TIME_Half(:,1:751);
            
            % Create Artificial Burst and Constructed Waveform for plotting
            ArtificialBurst = TIME_Half;
            
            ModelTitle = 'Half Cosine';
            
            clear half_ti half_burst SCALE_Half AddBaseline TIME_Half
                
        elseif strcmp(model(Modeli),'Gauss')
            
            % Get time
            Gauss_ti = -5:.001:5;
            
            % Fit Freq
            s=2./(2*pi.*BaseFrex);
            Gauss=exp(-Gauss_ti.^2./(2*s^2));
            Gauss=Gauss(1:4:end); % downsample
            
            % Amplitude scaling
            SCALE_Gauss=Gauss.*BaseScal;
            SCALE_Gauss=SCALE_Gauss(1,SCALE_Gauss(1,:)>0.001); % round tails for movement
            
            % Add Baseline and remove tail if beyond time window
            AddBaseline  = 251+BaseTime-111;
            TIME_Gauss = [zeros(1,AddBaseline),SCALE_Gauss,zeros(1,(751-(length(SCALE_Gauss)+AddBaseline)))];
            TIME_Gauss = TIME_Gauss(:,1:751);
            
            % Create Artificial Burst and Constructed Waveform for plotting
            ArtificialBurst = TIME_Gauss;
            
            ModelTitle = 'Gaussian';
            
            clear Gauss_ti s Gauss half_burst SCALE_Gauss AddBaseline TIME_Gauss
                    
        elseif strcmp(model(Modeli),'Wavelet')
            
            % Get time
            Wavelet_ti=-1:.001:1;
            
            % Fit Freq
            Delta=exp(2*1i*pi*BaseFrex.*Wavelet_ti);
            s=2./(2*pi.*BaseFrex);
            Delta_Taper=exp(-Wavelet_ti.^2./(2*s^2));
            
            % Fit Taper over the sine
            Delta_Taper_Shift=[zeros(1,BaseTapr),Delta_Taper(1:end-BaseTapr)];
            Deltalet=Delta.*Delta_Taper_Shift;
            
            % Manage temporal duration
            Deltalet=real(Deltalet(500:1500));    % 1001 samples
            % Re-sample based on srate of data
            Deltalet=resample(Deltalet,1,1000/srate);
            
            % Amplitude scaling
            Deltalet_scale=Deltalet.*BaseScal;
            Deltalet_scale=round(Deltalet_scale,3); % cut off tails to make fit, technically imprecise
            Deltalet_scale(:,Deltalet_scale(:,:)==0)=[]; % cut off tails to make fit, technically imprecise
            
            % Add Baseline and remove tail if beyond time window
            AddBaseline  = 61+BaseTime;
            TIME_Wavelet = [zeros(1,AddBaseline),Deltalet_scale,zeros(1,(751-(length(Deltalet_scale)+AddBaseline)))];
            TIME_Wavelet = TIME_Wavelet(:,1:751);
            
            % Create Artificial Burst and Constructed Waveform for plotting
            ArtificialBurst = TIME_Wavelet;
            
            ModelTitle = 'Wavelet';
            
            clear Wavelet_ti Delta s Delta_Taper Delta_Taper_Shift Deltalet Deltalet_scale AddBaseline TIME_Wavelet
                        
        elseif strcmp(model(Modeli),'Gamma')
            
            % Get time
            Gamma_ti = 0:.001:100;
            
            % Fit shape and rate params
            Gamma = makedist('Gamma','a',BaseShpe,'b',BaseRate);
            Gamma_Dist=pdf(Gamma,Gamma_ti);
            Gamma_Dist=Gamma_Dist(1:10:end); % downsample to fit with ERPs
            Gamma_Dist=Gamma_Dist(1,Gamma_Dist(1,:)>0.001);
            
            % Amplitude scaling
            SCALE_Gamma=Gamma_Dist.*BaseScal;
            
            % Add Baseline and remove tail if beyond time window
            AddBaseline  = 251+BaseTime+11;
            TIME_Gamma = [zeros(1,AddBaseline),SCALE_Gamma,zeros(1,(751-(length(SCALE_Gamma)+AddBaseline)))];
            TIME_Gamma = TIME_Gamma(:,1:751);
            
            % Create Artificial Burst and Constructed Waveform for plotting
            ArtificialBurst = TIME_Gamma;
            ModelTitle = 'Gamma';
            clear Gamma_ti Gamma Gamma_Dist SCALE_Gamma AddBaseline TIME_Gamma
            
        end
        
        % Artificial Waveforms
        subplot(size(model,2),3,PlotSpace(1:3)); hold on;
        plot(tx2disp,squeeze(ArtificialBurst),'c','linewidth',2);
        set(gca,'ylim',ylimit)
        title(ModelTitle);
        
        clear ArtificialBurst ModelTitle
        
end


%% Check GRID vs GRAD
for Modaliti = 1:2;
    load([homedir,Study,DataTitle,'HalfCosine_ModelFits_Grid']);
    GRID_params(Modaliti,:) = [mean([GRID_HALF(:,Modaliti).FREQ],2) mean([GRID_HALF(:,Modaliti).SCALING],2) round(mean([GRID_HALF(:,Modaliti).TIMESHIFT],2))];
    GRID_AIC(Modaliti,1) = mean([GRID_HALF(:,Modaliti).AIC],2);
    % fprintf([DataTitle,'\n GRID \n AIC = ', num2str(GRID_AIC),'\n Frequency = ',num2str(GRID_params(Modaliti,1),'\n
    load([homedir,Study,DataTitle,'HalfCosine_ModelFits_Grad']);
    GRAD_params(Modaliti,:) = [mean([GRAD_HALF(:,Modaliti).FREQ],2) mean([GRAD_HALF(:,Modaliti).SCALING],2) round(mean([GRAD_HALF(:,Modaliti).TIMESHIFT],2))];
    GRAD_AIC(Modaliti,1) = mean([GRAD_HALF(:,Modaliti).AIC],2);
end

%%
% Will put in correlations later
% 
%     for n = 1:size(Use.GRID_HALF,2)
%         Model_Half(n,1) = Use.GRID_HALF(n).FREQ;
%         Model_Half(n,2) = Use.GRID_HALF(n).SCALING;
%         Model_Half(n,3) = Use.GRID_HALF(n).TIMESHIFT;
% %         Model_Half(n,4) = Use.Age(n);
%         Model_Half(n,5) = Use.BDI(n);
%         Model_Half(n,6) = Use.BIS(n);
%         Model_Half(n,7) = Use.BAS_D(n);
%         Model_Half(n,8) = Use.BAS_FS(n);
%         Model_Half(n,9) = Use.BAS_RR(n);
%         Model_Half(n,10) = Use.MASQ_AD(n);
%         Model_Half(n,11) = Use.MASQ_GDD(n);
%         
%         Model_Full(n,1) = Use.GRID_FULL(n).FREQ;
%         Model_Full(n,2) = Use.GRID_FULL(n).SCALING;
%         Model_Full(n,3) = Use.GRID_FULL(n).TIMESHIFT;
% %         Model_Full(n,4) = Use.Age(n);
%         Model_Full(n,5) = Use.BDI(n);
%         Model_Full(n,6) = Use.BIS(n);
%         Model_Full(n,7) = Use.BAS_D(n);
%         Model_Full(n,8) = Use.BAS_FS(n);
%         Model_Full(n,9) = Use.BAS_RR(n);
%         Model_Full(n,10) = Use.MASQ_AD(n);
%         Model_Full(n,11) = Use.MASQ_GDD(n);
%         
%         Model_P3(n,1) = Use.GRID_P3(n).FREQ;
%         Model_P3(n,2) = Use.GRID_P3(n).SCALING;
%         Model_P3(n,3) = Use.GRID_P3(n).TIMESHIFT;
% %         Model_P3(n,4) = Use.Age(n);
%         Model_P3(n,5) = Use.BDI(n);
%         Model_P3(n,6) = Use.BIS(n);
%         Model_P3(n,7) = Use.BAS_D(n);
%         Model_P3(n,8) = Use.BAS_FS(n);
%         Model_P3(n,9) = Use.BAS_RR(n);
%         Model_P3(n,10) = Use.MASQ_AD(n);
%         Model_P3(n,11) = Use.MASQ_GDD(n);
%     end
%     
%     sizevars = 11;
%     LABELS = {'FREQ','SCALING','TIMESHIFT','AGE','BDI','BIS','BAS_D','BAS_FS','BAS_RR','MASQ_AD','MASQ_GDD'};
%     
% %     scatter(Use.BDI,Y)
%     
%     % Correlation Matrix for Half Cosine
%     [MATRIX_Half,Half_P]=corr(Model_Half,Model_Half,'rows','complete','Type','Spearman');
%     Half_P(Half_P<=.05)=NaN; Half_P(Half_P>.05)=0;  Half_P(isnan(Half_P))=1;
%     
%     figure; hold on;
%     imagesc(MATRIX_Half); axis ij;
%     set(gca,'clim',[-1 1],'xlim',[0 sizevars+1],'ylim',[0 sizevars+1],'xtick',[1:1:sizevars],'ytick',[1:1:sizevars],...
%         'xticklabels',LABELS,'yticklabels',LABELS);  
%     title('Model Parameters for Half Cosine');
%     cbar
%     
%     % Correlation Matrix for Full Cosine
%     [MATRIX_Full,Full_P]=corr(Model_Full,Model_Full,'rows','complete','Type','Spearman');
%     Full_P(Full_P<=.05)=NaN; Full_P(Full_P>.05)=0;  Full_P(isnan(Full_P))=1;
%     
%     figure; hold on;
%     imagesc(MATRIX_Full); axis ij;
%     set(gca,'clim',[-1 1],'xlim',[0 sizevars+1],'ylim',[0 sizevars+1],'xtick',[1:1:sizevars],'ytick',[1:1:sizevars],...
%         'xticklabels',LABELS,'yticklabels',LABELS);  
%     title('Model Parameters for Full Cosine');
%     cbar
%     
%     % Correlation Matrix for P3
%     [MATRIX_P3,P3_P]=corr(Model_P3,Model_P3,'rows','complete','Type','Spearman');
%     P3_P(P3_P<=.05)=NaN; P3_P(P3_P>.05)=0;  P3_P(isnan(P3_P))=1;
%     
%     figure; hold on;
%     imagesc(MATRIX_P3); axis ij;
%     set(gca,'clim',[-1 1],'xlim',[0 sizevars+1],'ylim',[0 sizevars+1],'xtick',[1:1:sizevars],'ytick',[1:1:sizevars],...
%         'xticklabels',LABELS,'yticklabels',LABELS);  
%     title('Model Parameters for P3 Half Cosine');
%     cbar
    
%     clearvars -except Iterations Modaliti homedir Training_Sight Training_Sound tx2disp Electrode srate FREQ_VECTOR TAPERSHIFT_VECTOR SCALING_VECTOR TIMESHIFT_VECTOR

    
% end

%% Set Training parameters aside

clearvars -except homedir Training_Sight Training_Sound tx2disp Electrode srate FREQ_VECTOR TAPERSHIFT_VECTOR SCALING_VECTOR TIMESHIFT_VECTOR

%% Fit to Other Data Sets
% Outputs must be: [Pun_ERP    Rew_ERP    RewP_Diff    tx2disp   srate   DataTitle]
% Where each ERP is a 1D vector, tx2disp is the samples-to-time conversion, srate is sampling rate, DataTitle is for saving the output

Datasets = 6;

for Dataseti = 6%2:Datasets
    
    
    %% Load the data
    
%     if Dataseti == 1 % Single Modality Reward
%         Iterations = 2;
%         rootdir = [homedir,'2_SightOrSound\'];
%         cd(rootdir);
%         
%         % Load Sight FCz
%         load('Load_Sight_SingleModalityReward.mat');
%         Sight.Pun_ERP = Pun_ERP; Sight.Rew_ERP = Rew_ERP;
%         Sight.RewP_Diff = RewP_Diff; Sight.DataTitle = DataTitle;
%         Sight.GraphTitle = 'Exp1 - SoS Sight';
%         clear Pun_ERP Rew_ERP RewP_Diff DataTitle
%         
%         % Load Sound FCz
%         load('Load_Sound_SingleModalityReward.mat');
%         Sound.Pun_ERP = Pun_ERP; Sound.Rew_ERP = Rew_ERP;
%         Sound.RewP_Diff = RewP_Diff; Sound.DataTitle = DataTitle;
%         Sound.GraphTitle = 'Exp1 - SoS Sound';
%         clear Pun_ERP Rew_ERP RewP_Diff DataTitle
        
    if Dataseti == 2 % Single Modality Reward with P3 (Darin Study 1)
        Iterations = 2;
        rootdir = [homedir,'3_SightOrSoundP3\'];
        cd(rootdir);
        
        % Load Data
        load('V_SoS_30-Aug-2019.mat');
        
        % Pull Sight Data
        Sight.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,4),1) ); 
        Sight.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2),1) );
        Sight.RewP_Diff = Sight.Rew_ERP - Sight.Pun_ERP; 
        Sight.DataTitle = 'Sight_SingleModalityReward_P3';
        Sight.GraphTitle = 'Exp2 - SoS_P3 Sight';
        Sight.time1 = 200;
        Sight.time2 = 500;
        
        % Pull Sound Data
        Sound.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1),1) ); 
        Sound.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3),1) );
        Sound.RewP_Diff = Sound.Rew_ERP - Sound.Pun_ERP; 
        Sound.DataTitle = 'Sound_SingleModalityReward_P3';
        Sound.GraphTitle = 'Exp2 - SoS_P3 Sound';
        Sound.time1 = 0;
        Sound.time2 = 300;
        
        tx2disp=-500:2:1000;
        clear MEGA_ERPs ROI
    
    elseif Dataseti == 3 % Dual Modality Reward
        Iterations = 3;
        rootdir = [homedir,'5_SightAndSound\'];
        cd(rootdir);
        
        % Load Sight FCz
        load('Load_Sight_DualModalityReward.mat');
        Sight.Pun_ERP = Pun_ERP; Sight.Rew_ERP = Rew_ERP;
        Sight.RewP_Diff = RewP_Diff; Sight.DataTitle = DataTitle;
        Sight.GraphTitle = 'Exp3 - SnS Sight';
        Sight.time1 = 200;
        Sight.time2 = 500;
        clear Pun_ERP Rew_ERP RewP_Diff DataTitle
        
        % Load Sound FCz
        load('Load_Sound_DualModalityReward.mat');
        Sound.Pun_ERP = Pun_ERP; Sound.Rew_ERP = Rew_ERP;
        Sound.RewP_Diff = RewP_Diff; Sound.DataTitle = DataTitle;
        Sound.GraphTitle = 'Exp3 - SnS Sound';
        Sound.time1 = 0;
        Sound.time2 = 300;
        clear Pun_ERP Rew_ERP RewP_Diff DataTitle
        
        % Load Combo FCz
        load('Load_Combo_DualModalityReward.mat');
        Combo.Pun_ERP = Pun_ERP; Combo.Rew_ERP = Rew_ERP;
        Combo.RewP_Diff = RewP_Diff; Combo.DataTitle = DataTitle;
        Combo.GraphTitle = 'Exp3 - SnS Combo';
        Combo.time1 = 0;
        Combo.time2 = 300;
        Combo.time3 = 200;
        Combo.time4 = 500;
        clear Pun_ERP Rew_ERP RewP_Diff DataTitle
        
    elseif Dataseti == 4 % Dual Modality Reward with P3 (Darin Study 2)
        Iterations = 3;
        rootdir = [homedir,'6_SightAndSoundP3\'];
        cd(rootdir);
        
        load('V_SnS_30-Aug-2019.mat');
        
        % Pull Sight Data
        Sight.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1),1) ); 
        Sight.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,3),1) );
        Sight.RewP_Diff = Sight.Rew_ERP - Sight.Pun_ERP; 
        Sight.DataTitle = 'Sight_DualModalityReward_P3';
        Sight.GraphTitle = 'Exp4 - SnS_P3 Sight';
        Sight.time1 = 200;
        Sight.time2 = 500;
        
        % Pull Sound Data
        Sound.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1),1) ); 
        Sound.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,4),1) );
        Sound.RewP_Diff = Sound.Rew_ERP - Sound.Pun_ERP; 
        Sound.DataTitle = 'Sound_DualModalityReward_P3';
        Sound.GraphTitle = 'Exp4 - SnS_P3 Sound';
        Sound.time1 = 0;
        Sound.time2 = 300;
        
        % Pull Combo Data
        Combo.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1),1) ); 
        Combo.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2),1) );
        Combo.RewP_Diff = Combo.Rew_ERP - Combo.Pun_ERP; 
        Combo.DataTitle = 'Combo_DualModalityReward_P3';
        Combo.GraphTitle = 'Exp4 - SnS_P3 Combo';
        Combo.time1 = 0;
        Combo.time2 = 300;
        
        tx2disp=-500:2:1000;
        clear MEGA_ERPs ROI
        
    elseif Dataseti == 5 % Sound Delayed Reward (GREMLIN)
        YLim = [-14 4]; 
        Iterations = 3;
        rootdir = [homedir,'7_SoundDelay_Gremlin\'];
        cd(rootdir);
        
        load('V_Gremlin_07-Sep-2021.mat');
        % Pull 100 Data
        Shrt.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1,1),1) ); 
        Shrt.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,1),1) );
        Shrt.RewP_Diff = Shrt.Rew_ERP - Shrt.Pun_ERP; 
        Shrt.time1 = 100; Shrt.time2 = 600;
        Shrt.DataTitle = 'Shrt_SoundDelay';
        Shrt.GraphTitle = '100ms Delay';
        
        % Pull 300 Data
        Medi.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1,2),1) ); 
        Medi.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,2),1) );
        Medi.RewP_Diff = Medi.Rew_ERP - Medi.Pun_ERP;
        Medi.time1 = 400; Medi.time2 = 800;
        Medi.DataTitle = 'Medi_SoundDelay';
        Medi.GraphTitle = '300ms Delay';
        Medi.TIMESHIFT_VECTOR = 100:1:300;
        
        % Pull 500 Data
        Long.Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,1,3),1) ); 
        Long.Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,Electrode,1:751,2,3),1) );
        Long.RewP_Diff = Long.Rew_ERP - Long.Pun_ERP;
        Long.time1 = 600; Long.time2 = 1000;
        Long.DataTitle = 'Long_SoundDelay';
        Long.GraphTitle = '500ms Delay';
        Long.TIMESHIFT_VECTOR = 200:1:400;
        
        tx2disp=-500:2:1000;
        clear MEGA_ERPs USEsubjs
        
    elseif Dataseti == 6 % Visual Door Depressed
        Iterations = 1;
        rootdir = [homedir,'8_SightDep_R01Door\'];
        cd(rootdir);
        
        % ##### Soft Coded Stuff ##### %
        Electrode = 36;
        % ############################ %
        
        % Load Sight
        load('R01DepDoorVolt_28-Apr-2022.mat');
        
        % Individual differences
        Sight.Pun_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,:,3) );
        Sight.Rew_ERP = squeeze( DoorVolt.ERPs_FB(:,Electrode,:,1) );
        
        % Everything Else
        Sight.RewP_Diff = Sight.Rew_ERP - Sight.Pun_ERP;
%         Sight.Age = DoorVolt.Age;
        Sight.BDI = DoorVolt.BDI;
        Sight.BIS = DoorVolt.BIS;
        Sight.BAS_D = DoorVolt.BAS_D;
        Sight.BAS_FS = DoorVolt.BAS_FS;
        Sight.BAS_RR = DoorVolt.BAS_RR;
        Sight.MASQ_AD = DoorVolt.MASQ_AD;
        Sight.MASQ_GDD = DoorVolt.MASQ_GDD;
        
        srate = 500;
        DataTitle = 'R01_DepressedDoor_';
        Iterations = 1;
        Sight.DataTitle = DataTitle;
        Sight.GraphTitle = 'Exp6-Dep Doors';
        Sight.time1 = 200;
        Sight.time2 = 800;
        clear Pun_ERP Rew_ERP RewP_Diff DataTitle
        
    end
    

    %% Fit parameters to Training with scale and time as open parameters.
    
    for Iteri = 1:Iterations
        if Iteri == 1 && Dataseti ~= 5
            Use = Sight;
            Training = Training_Sight;
            ModalityTitle = 'Sight';
        elseif Iteri == 2 && Dataseti ~= 5
            Use = Sound;
            Training = Training_Sound;
            ModalityTitle = 'Sound';
        elseif Iteri == 3 && Dataseti ~= 5
            Use = Combo;
            Training = Training_Sound;
            Training2 = Training_Sight;
            ModalityTitle = 'Combo';
        elseif Iteri == 1 && Dataseti == 5
            Use = Shrt;
            TIMESHIFT_VECTOR = 0:1:200;
            Training = Training_Sound;
            ModalityTitle = 'Short';
        elseif Iteri == 2 && Dataseti == 5
            Use = Medi;
            TIMESHIFT_VECTOR = Medi.TIMESHIFT_VECTOR;
            Training = Training_Sound;
            ModalityTitle = 'Medium';
        elseif Iteri == 3 && Dataseti == 5
            Use = Long;
            TIMESHIFT_VECTOR = Long.TIMESHIFT_VECTOR;
            Training = Training_Sound;
            ModalityTitle = 'Long';
        end
        
        tic
        % Get Frequencies for separate models and compile for Freq Imput
        % here
        for n = 1:size(Use.RewP_Diff,2) %" if running by participant. Add it back for single. 
        
        FrequenciesFound = [find(FREQ_VECTOR==Training.GRID_HALF(n).FREQ) find(FREQ_VECTOR==Training.GRID_FULL(n).FREQ) find(FREQ_VECTOR==Training.GRID_P3(n).FREQ)];
        
        for FREQi=1:length(FrequenciesFound)
            for SCALINGi=1:length(SCALING_VECTOR)
                for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                    
                    params=[FREQ_VECTOR( FrequenciesFound(FREQi) ),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR(TIMESHIFTi)];
                    
%                     [SSE_Half(FREQi,SCALINGi,TIMESHIFTi) SSE_Full(FREQi,SCALINGi,TIMESHIFTi) SSE_P3(FREQi,SCALINGi,TIMESHIFTi)] = MakeRewP_P3(params,Use.RewP_Diff,tx2disp,srate,Use.time1,Use.time2);
                    [SSE_Half(FREQi,SCALINGi,TIMESHIFTi),SSE_Full(FREQi,SCALINGi,TIMESHIFTi),SSE_P3(FREQi,SCALINGi,TIMESHIFTi)] = MakeRewP_P3(params,Use.RewP_Diff(n,:),tx2disp,srate,Use.time1,Use.time2);    
                    
                    clear params;
                    
                end
            end
        end
        
        if Iteri == 3 && Dataseti ~= 5
            FrequenciesFound = [find(FREQ_VECTOR==Training.GRID_HALF.FREQ) find(FREQ_VECTOR==Training.GRID_FULL.FREQ) find(FREQ_VECTOR==Training.GRID_P3.FREQ)];
            
            for FREQi=1:length(FrequenciesFound)
                for SCALINGi=1:length(SCALING_VECTOR)
                    for TIMESHIFTi=1:length(TIMESHIFT_VECTOR)
                        
                        params=[FREQ_VECTOR( FrequenciesFound(FREQi) ),SCALING_VECTOR(SCALINGi),TIMESHIFT_VECTOR(TIMESHIFTi)];
                        
                        [SSE_Half(FREQi,SCALINGi,TIMESHIFTi) SSE_Full(FREQi,SCALINGi,TIMESHIFTi) SSE_P3(FREQi,SCALINGi,TIMESHIFTi)] = MakeRewP_P3(params,Use.RewP_Diff,tx2disp,srate,Use.time1,Use.time2);
                        
                        
                        clear params;
                        
                    end
                end
            end
        end
            
            
        
        GRID_HALF.toc=toc;
        
        Use_SSE_Half = squeeze(SSE_Half(1,:,:));
        Use_SSE_Full = squeeze(SSE_Full(2,:,:));
        Use_SSE_P3 = squeeze(SSE_P3(3,:,:));
        
        
        % ------------------------------------- find the values that correspond
        % to the lowest SSE for a half cosine
        
        dims=size(Use_SSE_Half);
        [lowest_error,lowest_idx]=min(Use_SSE_Half(Use_SSE_Half(:)>0));
        [B,C] = ind2sub(dims, lowest_idx);
        
        % ########
        GRID_HALF.SSE=lowest_error;
        GRID_HALF.FREQ=Training.GRID_HALF.FREQ;
        GRID_HALF.SCALING=SCALING_VECTOR(B);
        GRID_HALF.TIMESHIFT=TIMESHIFT_VECTOR(C);
        % ########
        
        % ------------------------------------- find the values that correspond
        % to the lowest SSE for a full cosine
        dims=size(Use_SSE_Full);
        [lowest_error,lowest_idx]=min(Use_SSE_Full(Use_SSE_Full(:)>0));
        [B,C] = ind2sub(dims, lowest_idx);
        
        % ########
        GRID_FULL.SSE=lowest_error;
        GRID_FULL.FREQ=Training.GRID_FULL.FREQ;
        GRID_FULL.SCALING=SCALING_VECTOR(B);
        GRID_FULL.TIMESHIFT=TIMESHIFT_VECTOR(C);
        % ########
        
        % ------------------------------------- find the values that correspond
        % to the lowest SSE for a full cosine
        dims=size(Use_SSE_P3);
        [lowest_error,lowest_idx]=min(Use_SSE_P3(Use_SSE_P3(:)>0));
        [B,C] = ind2sub(dims, lowest_idx);
        
        % ########
        GRID_P3.SSE=lowest_error;
        GRID_P3.FREQ=Training.GRID_P3.FREQ;
        GRID_P3.SCALING=SCALING_VECTOR(B);
        GRID_P3.TIMESHIFT=TIMESHIFT_VECTOR(C);
        % ########
    

        Testing.SSE_Half(n,:,:,:) = GRID_HALF;
        Testing.SSE_Full(n,:,:,:) = GRID_FULL;
        Testing.SSE_P3(n,:,:,:) = GRID_P3;
        
%         %% Create the best-fitting RewP difference wave and Save Images
%                 
%         % ########
%         half_FREQ_fit= GRID_HALF.FREQ;
%         half_SCALE_fit=GRID_HALF.SCALING;
%         half_TSHIFT_fit=round(GRID_HALF.TIMESHIFT);
%         
%         full_FREQ_fit= GRID_FULL.FREQ;
%         full_SCALE_fit=GRID_FULL.SCALING;
%         full_TSHIFT_fit=round(GRID_FULL.TIMESHIFT);
%         
%         P3_FREQ_fit= GRID_P3.FREQ;
%         P3_SCALE_fit=GRID_P3.SCALING;
%         P3_TSHIFT_fit=round(GRID_P3.TIMESHIFT);
%         
%         % ########
%         
%         timedur=length(tx2disp);
%         
%         % % ######## new code to compare full vs. half cosine wave
%         half_ti = -((1/half_FREQ_fit)*.25):(1000/srate)/1000:((1/half_FREQ_fit)*.25);
%         full_ti = -((1/full_FREQ_fit)*.25):(1000/srate)/1000:((1/full_FREQ_fit)*.75);
%         P3_ti   = ((1/P3_FREQ_fit)*.25):(1000/srate)/1000:((1/P3_FREQ_fit)*.75);
%         
%         % Fit Freq
%         half_burst = cos(2*pi*half_FREQ_fit.*half_ti);
%         full_burst = cos(2*pi*full_FREQ_fit.*full_ti);
%         P3_burst   = cos(2*pi*P3_FREQ_fit.*P3_ti);
%         
%         % Amplitude scaling
%         SCALE_Half = half_burst.*half_SCALE_fit;
%         SCALE_Full = full_burst.*full_SCALE_fit;
%         SCALE_P3   = P3_burst.*P3_SCALE_fit;
%         
%         % Shift the RewP by time
%         AddBaseline_half = 251+half_TSHIFT_fit;
%         AddBaseline_full = 251+full_TSHIFT_fit;
%         AddBaseline_P3 = 251+P3_TSHIFT_fit+100;
%         SimRewP_Half=[zeros(1,AddBaseline_half),SCALE_Half,zeros(1,(751-(length(SCALE_Half)+AddBaseline_half)))];
%         SimRewP_Half=SimRewP_Half(:,1:751);
%         SimRewP_Full=[zeros(1,AddBaseline_full),SCALE_Full,zeros(1,(751-(length(SCALE_Full)+AddBaseline_full)))];
%         SimRewP_Full=SimRewP_Full(1,1:751);
%         SimRewP_P3   = [zeros(1,AddBaseline_P3), SCALE_P3 ,zeros(1,(751-(length( SCALE_P3 )+AddBaseline_P3)))];
%         SimRewP_P3   = SimRewP_P3(1,1:751);
%         % ########
%         T1 = find(tx2disp(1,:)==Use.time1); T2 = find(tx2disp(1,:)==Use.time2);
%         
%         %     SSE_Combined_Half_P3 = sum((Use.RewP_Diff(n,T1:T2)-( SimRewP_Half(1,T1:T2) + SimRewP_P3(1,T1:T2))).^2);
%         SSE_Combined_Half_P3 = sum((Use.RewP_Diff(T1:T2)'-( SimRewP_Half(1,T1:T2) + SimRewP_P3(1,T1:T2))).^2);
%         
%         figure;  hold on;
%         subplot(3,3,1); hold on;
%         text(.1,.8,['SSE: ',num2str(GRID_HALF.SSE)],'sc');
%         text(.1,.65,['FREQ: ',num2str(GRID_HALF.FREQ)],'sc');
%         text(.1,.5,['SCAL: ',num2str(GRID_HALF.SCALING)],'sc');
%         text(.1,.35,['TIME: ',num2str(GRID_HALF.TIMESHIFT)],'sc');
%         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
%         subplot(3,3,4); hold on;
%         text(.1,.8,['SSE: ',num2str(GRID_FULL.SSE)],'sc');
%         text(.1,.65,['FREQ: ',num2str(GRID_FULL.FREQ)],'sc');
%         text(.1,.5,['SCAL: ',num2str(GRID_FULL.SCALING)],'sc');
%         text(.1,.35,['TIME: ',num2str(GRID_FULL.TIMESHIFT)],'sc');
%         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
%         subplot(3,3,7); hold on;
%         text(.1,.8,['SSE: ',num2str(SSE_Combined_Half_P3)],'sc');
%         text(.1,.65,['FREQ: ',num2str(GRID_P3.FREQ)],'sc');
%         text(.1,.5,['SCAL: ',num2str(GRID_P3.SCALING)],'sc');
%         text(.1,.35,['TIME: ',num2str(GRID_P3.TIMESHIFT)],'sc');
%         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
%         
% %         text(.1,.8,['SSE: ',num2str(GRID_HALF(n).SSE)],'sc');
% %         text(.1,.65,['FREQ: ',num2str(GRID_HALF(n).FREQ)],'sc');
% %         text(.1,.5,['SCAL: ',num2str(GRID_HALF(n).SCALING)],'sc');
% %         text(.1,.35,['TIME: ',num2str(GRID_HALF(n).TIMESHIFT)],'sc');
% %         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
% %         subplot(3,3,4); hold on;
% %         text(.1,.8,['SSE: ',num2str(GRID_FULL(n).SSE)],'sc');
% %         text(.1,.65,['FREQ: ',num2str(GRID_FULL(n).FREQ)],'sc');
% %         text(.1,.5,['SCAL: ',num2str(GRID_FULL(n).SCALING)],'sc');
% %         text(.1,.35,['TIME: ',num2str(GRID_FULL(n).TIMESHIFT)],'sc');
% %         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
% %         subplot(3,3,7); hold on;
% %         text(.1,.8,['SSE: ',num2str(SSE_Combined_Half_P3)],'sc');
% %         text(.1,.65,['FREQ: ',num2str(GRID_P3(n).FREQ)],'sc');
% %         text(.1,.5,['SCAL: ',num2str(GRID_P3(n).SCALING)],'sc');
% %         text(.1,.35,['TIME: ',num2str(GRID_P3(n).TIMESHIFT)],'sc');
% %         text(.1,.2,['DATA: ',Use.GraphTitle],'sc');
%         
%         subplot(3,3,2); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Pun','Rew','Diff'},'location','NorthWest')
%         title([ModalityTitle,' Real ERPs']);
%         subplot(3,3,5); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Pun','Rew','Diff'},'location','NorthWest')
%         title([ModalityTitle,' Real ERPs']);
%         subplot(3,3,8); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Pun','Rew','Diff'},'location','NorthWest')
%         title([ModalityTitle,' Real ERPs']);
%         
%         subplot(3,3,3); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,SimRewP_Half,'c','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_Half,'k','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,SimRewP_Half,'c','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_Half,'k','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Pun','SimRewP','Combo'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP Half']);
%         subplot(3,3,6); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,SimRewP_Full,'c','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_Full,'k','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,SimRewP_Full,'c','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_Full,'k','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         
%         set(gca,'ylim',YLim)
%         legend({'Pun','SimRewP','Combo'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP Full']);
%         subplot(3,3,9); hold on;
%         %     plot(tx2disp,Use.Pun_ERP(n,:),'r','linewidth',2);
%         %     plot(tx2disp,SimRewP_P3,'c','linewidth',2);
%         %     plot(tx2disp,SimRewP_Half,'m','linewidth',2);
%         % %     plot(tx2disp,Use.Pun_ERP'+SimRewP_P3,'k','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_P3+SimRewP_Half,'k-','linewidth',2);
%         %     plot(tx2disp,Use.RewP_Diff(n,:),'b','linewidth',2);
%         
%         plot(tx2disp,Use.Pun_ERP','r','linewidth',2);
%         plot(tx2disp,SimRewP_P3,'c','linewidth',2);
%         plot(tx2disp,SimRewP_Half,'m','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP'+SimRewP_P3,'k','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_P3+SimRewP_Half,'k-','linewidth',2);
%         plot(tx2disp,Use.RewP_Diff','b','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Pun','SimRewP','Combo'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP P3 + Half']);
%         
%         figure;
%         subplot(1,3,1); hold on;
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,SimRewP_Half,'c','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_Half,'k','linewidth',2);
%         
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,SimRewP_Half,'c','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_Half,'k','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Rew','SimRewP','Combo - Pun+Sim'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP vs. Real RewP']);
%         subplot(1,3,2); hold on;
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,SimRewP_Full,'c','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_Full,'k','linewidth',2);
%         
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,SimRewP_Full,'c','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_Full,'k','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Rew','SimRewP','Combo - Pun+Sim'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP vs. Real RewP']);
%         subplot(1,3,3); hold on;
%         %     plot(tx2disp,Use.Rew_ERP(n,:),'g','linewidth',2);
%         %     plot(tx2disp,SimRewP_P3,'c','linewidth',2);
%         %     plot(tx2disp,SimRewP_Half,'m','linewidth',2);
%         %     plot(tx2disp,Use.Pun_ERP(n,:)+SimRewP_P3+SimRewP_Half,'k-','linewidth',2);
%         
%         plot(tx2disp,Use.Rew_ERP','g','linewidth',2);
%         plot(tx2disp,SimRewP_P3,'c','linewidth',2);
%         plot(tx2disp,SimRewP_Half,'m','linewidth',2);
%         plot(tx2disp,Use.Pun_ERP'+SimRewP_P3+SimRewP_Half,'k','linewidth',2);
%         set(gca,'ylim',YLim)
%         legend({'Rew','SimRewP','Combo - Pun+Sim'},'location','NorthWest')
%         title([ModalityTitle,' Real Pun + SimRewP vs. Real RewP']);
%         
%         Experiment{Dataseti}.Iteration{Iteri,1} = GRID_HALF;
%         Experiment{Dataseti}.Iteration{Iteri,2} = GRID_FULL;
% %         Experiment{Dataseti}.Iteration{Iteri,3} = SSE_Half_Match_Full;
% %         Experiment{Dataseti}.Iteration{Iteri,4} = SSE_Full_Match_Half;
%         Experiment{Dataseti}.Iteration{Iteri,5} = Use;
        

        end
        
        clear Use Training
        
    end
    
    clearvars -except Experiment homedir Training_Sight Training_Sound tx2disp Electrode srate time1 time2 FREQ_VECTOR TAPERSHIFT_VECTOR SCALING_VECTOR TIMESHIFT_VECTOR

end


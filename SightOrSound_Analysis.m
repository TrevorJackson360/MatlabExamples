%% Sight and Sound Analysis

close all; clear all; clc

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Set Up Directories and BIDS Stuff
rootdir=('G:\UNM_Data\SightOrSound2\');
sourcedir=[rootdir,'sourcedata\'];
scriptdir = [rootdir,'Scripts\'];
locpath=('Y:\Programs\eeglab12_0_2_1b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp');

% addpath(['Y:\Programs\Bidsify\']);
% addpath(genpath('Y:\Programs\eeglab14_1_2b'));
addpath(scriptdir);

% ################
TASKNAME = 'SOS2';
FileFormat = 'BV';
DATE = date;
ElectrodesForTF = {'Fz','FCz','Cz','CPz'};
% ################

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get Participants and Remove Any Bad Ones (If Necessary)

cd(sourcedir);
sx = dir('*sub-*');
for si=1:length(sx)
    subjs(si)=str2num(sx(si).name(5:end));
end

% Get rid of Bad Participants for whatever reason
BadParticipants = [10059];
if ~isempty(BadParticipants)
    
    BadPartLoc = zeros(1,size(BadParticipants, 2));
    
    % Baseline
    for i = 1:size(BadParticipants,2)
        if ~isempty(find(subjs(1,:)==BadParticipants(1,i)))
            BadPartLoc(1,i) = find(subjs(1,:)==BadParticipants(1,i));
        end
    end
    clear i
    BadPartLoc = sort(BadPartLoc, 'descend');
    for i = size(BadPartLoc,2):-1:1
        if BadPartLoc(1,i)==0
            BadPartLoc(i)=[];
        end
    end
    BadPartLoc = sort(BadPartLoc, 'ascend');
    clear i
    for i = size(BadPartLoc,2):-1:1
        subjs(BadPartLoc(i)) = [];
        sx(BadPartLoc(i)) = [];
    end
    si=size(subjs,2);
    clear i BadPartLoc
end


%% Clean the Data

Step1_EEG_CLEANUP_SoS;

%% Save the Data for RewP Model
% Outputs must be: [Pun_ERP  Rew_ERP  RewP_Diff  tx2disp  srate  DataTitle]
% Where each ERP is a 1D vector, tx2disp is the samples-to-time conversion,
% srate is sampling rate, DataTitle is for saving the output.

load BV_ChanLocs_60
ChannelForModel = 'FCz';
chani = find(strcmpi({ChannelForModel},{BV_Chanlocs_60.labels}));

% Sight Specific
Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,chani,1:751,4),1) );
Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,chani,1:751,2),1) );
RewP_Diff = Rew_ERP-Pun_ERP;
tx2disp = -500:2:1000;
srate = 500;
DataTitle = 'Sight_SingleModalityReward';
save(['Y:\EEG_Data\JACKSON\RewP Model\ModalityReward\Load_',DataTitle],...
    'Pun_ERP','Rew_ERP','RewP_Diff','tx2disp','srate','DataTitle');
clear Pun_ERP Rew_ERP RewP_Diff DataTitle

% Sound Specific
Pun_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,chani,1:751,1),1) );
Rew_ERP = squeeze( mean(MEGA_ERPs.FBbyCUE(:,chani,1:751,3),1) );
RewP_Diff = Rew_ERP-Pun_ERP;
DataTitle = 'Sound_SingleModalityReward';
save(['Y:\EEG_Data\JACKSON\RewP Model\ModalityReward\Load_',DataTitle],...
    'Pun_ERP','Rew_ERP','RewP_Diff','tx2disp','srate','DataTitle');
clear Pun_ERP Rew_ERP RewP_Diff DataTitle


%% Plot the Data

Step2_Plot_TF_and_ERPs_SoS;

%% Get Weird with the Data

% Step4_PCA;

%% Pull Stuff for SPSS file

TrevorSaysYes = 0;

if TrevorSaysYes == 1;
    load('G:\UNM_Data\RewP Model\ModalityReward\E2_SightOrSound\SoS_Data.mat');
end
% ERP Stuff
tx2disp=-500:2:1000;
T1 = 200; T2 = 350;
T3 = 100; T4 = 250;
time1 = find(tx2disp==T1); time2 = find(tx2disp==T2);
time3 = find(tx2disp==T3); time4 = find(tx2disp==T4);
frex=logspace(.01,1.7,50);
frex1 = 1; frex2 = 18; frex3 = 19; frex4 = 27;
SiRw = 2; SiPn = 4; SoRw = 3; SoPn = 1;

% ERPs (FCz)
MEGA.Sight_FRN = [squeeze(mean(MEGA_ERPs.FBbyCUE(:,36,time1:time2,SiPn),3))];
MEGA.Sound_FRN = [squeeze(mean(MEGA_ERPs.FBbyCUE(:,36,time3:time4,SoPn),3))];
MEGA.Sight_RwP = [squeeze(mean(MEGA_ERPs.FBbyCUE(:,36,time1:time2,SiRw),3))];
MEGA.Sound_RwP = [squeeze(mean(MEGA_ERPs.FBbyCUE(:,36,time3:time4,SoRw),3))];
MEGA.Sight_Dif = [squeeze(mean(  (MEGA_ERPs.FBbyCUE(:,36,time1:time2,SiRw) - MEGA_ERPs.FBbyCUE(:,36,time1:time2,SiPn))  ,3))];
MEGA.Sound_Dif = [squeeze(mean(  (MEGA_ERPs.FBbyCUE(:,36,time3:time4,SoRw) - MEGA_ERPs.FBbyCUE(:,36,time3:time4,SoPn))  ,3))];

% TF (FCz)
MEGA.TFN_Si = squeeze( mean(  mean( MEGA_TF.FB(:,2,frex3:frex4,time1:time2,SiPn) ,3)  ,4));
MEGA.DRP_Si = squeeze( mean(  mean( MEGA_TF.FB(:,2,frex1:frex2,time1:time2,SiRw) ,3)  ,4));
MEGA.TFN_So = squeeze( mean(  mean( MEGA_TF.FB(:,2,frex3:frex4,time3:time4,SoPn) ,3)  ,4));
MEGA.DRP_So = squeeze( mean(  mean( MEGA_TF.FB(:,2,frex1:frex2,time3:time4,SoRw) ,3)  ,4));

MEGA.ITPC_TFRN_Si = squeeze( mean(  mean( MEGA_ITPC.FB(:,2,frex3:frex4,time1:time2,SiPn) ,3)  ,4));
MEGA.ITPC_DRew_Si = squeeze( mean(  mean( MEGA_ITPC.FB(:,2,frex1:frex2,time1:time2,SiRw) ,3)  ,4));
MEGA.ITPC_TFRN_So = squeeze( mean(  mean( MEGA_ITPC.FB(:,2,frex3:frex4,time3:time4,SoPn) ,3)  ,4));
MEGA.ITPC_DRew_So = squeeze( mean(  mean( MEGA_ITPC.FB(:,2,frex1:frex2,time3:time4,SoRw) ,3)  ,4));

SPSS(:,1) = MEGA.Sight_FRN;
SPSS(:,2) = MEGA.Sound_FRN;
SPSS(:,3) = MEGA.Sight_RwP;
SPSS(:,4) = MEGA.Sound_RwP;
SPSS(:,5) = MEGA.Sight_Dif;
SPSS(:,6) = MEGA.Sound_Dif;
SPSS(:,7) = MEGA.TFN_Si;
SPSS(:,8) = MEGA.DRP_Si;
SPSS(:,9) = MEGA.TFN_So;
SPSS(:,10) = MEGA.DRP_So;
SPSS(:,11) = MEGA.ITPC_TFRN_Si;
SPSS(:,12) = MEGA.ITPC_DRew_Si;
SPSS(:,13) = MEGA.ITPC_TFRN_So;
SPSS(:,14) = MEGA.ITPC_DRew_So;


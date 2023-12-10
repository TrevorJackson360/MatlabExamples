%% Step 1: Clean EEG / Set up ERP, TF, and ITPC

CB1 = 1;
CB2 = 1;
ALL = 1;

for subji=1:size(subjs,2)
    subno=num2str(subjs(subji));
    clc; disp(['TRYING  Subj: ',subno])
    cd(sourcedir);
    
    if  ~exist([rootdir,'derivatives\eegpreprocess\sub-',subno,'\eeg\',num2str(subno),'_',TASKNAME,'_ready.mat'])
        
        % Random soft-coding in case mastoids get mixed with eyes
        if subno == '10003'
            MastoidChannels = [10 27];
            EyeChannels = [5 21];
        else
            MastoidChannels = [10 21];
            EyeChannels = [5 27];
        end
        
        % Load Data
        clc;  format shortg
        disp(clock);
        
        eeglab;
        
        EEG = pop_loadbv([sourcedir,'sub-',subno,'\eeg\'],[subno,'_',TASKNAME,'.vhdr']);
        disp(['LOADED  Subj: ',subno])
        
        % Get Locs
        EEG = pop_chanedit(EEG,    'lookup', locpath);
        EEG = eeg_checkset( EEG );
        EEG.chanlocs(MastoidChannels(1)).labels = 'M1'; EEG.chanlocs(MastoidChannels(2)).labels = 'M2';
        % Get event types
        for ai=2:length(EEG.event); clear temp; temp=EEG.event(ai).type;
            if isempty(strmatch('boundary',temp)); TYPES(ai)=str2num(temp(2:end)) ; clear temp; end
        end
        UNIQUE_TYPES=unique(TYPES);
        for ai=1:length(UNIQUE_TYPES); UNIQUE_TYPES_COUNT(ai)=sum(TYPES==UNIQUE_TYPES(ai)); end
        clc; TRIGGERS=[UNIQUE_TYPES;UNIQUE_TYPES_COUNT] % Trigger type, Frequency
        clear UNIQUE_TYPES UNIQUE_TYPES_COUNT TRIGGERS
        
        
        
        
        % -----------  Triggers
        % 10 = Door Presentation
        % 15X = Button Press
        %       X = 1: Left Press
        %       X = 2: Middle Press
        %       X = 3: Right Press
        % 1 = Sight Reward
        % 2 = Sound Reward
        % 3 = Sight Null
        % 4 = Sound Null
        
        % Epoch
        All_STIM= {'S 10'};
        EEG = pop_epoch( EEG, All_STIM, [-2         6], 'newname', 'Epochs', 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        
        % EVENT COLUMNS!!!
        EVENTSNames = {'TrialNum','FB','FB_Latency','Simple_FB'};
        
        % EVENTS COLUMNS!!!
        % 1 - Trial Num
        % 2 - FB
        % 3 - FB Latency
        % 4 - FB as 1 (Reward) or 0 (Null)
        
        % Get the good info out of the epochs
        FB=[1,2,3,4];
        for ai=1:size(EEG.epoch,2)
            EEG.epoch(ai).FB=NaN; EEG.epoch(ai).FBlat=NaN;
            for bi=1:size(EEG.epoch(ai).eventlatency,2)
                if EEG.epoch(ai).eventlatency{bi}==0
                    if size(EEG.epoch(ai).eventtype,2)>2
                        tempname=str2num(EEG.epoch(ai).eventtype{bi+2}(2:end));
                        if any(tempname==FB)
                            clear tempname;
                            tempname=str2num(EEG.epoch(ai).eventtype{bi+2}(2:end));
                            EVENTS(ai,1)=ai;
                            EVENTS(ai,2)=tempname;
                            EVENTS(ai,3) = EEG.epoch(ai).eventlatency{bi+2};
                            if tempname == 1 || tempname == 2 
                                EVENTS(ai,4) = 1;
                            elseif tempname == 4 || tempname == 3
                                EVENTS(ai,4) = 0;
                            end
                            EEG.epoch(ai).FB=tempname;  clear tempname;
                        end
                    else
                        EVENTS(ai,1)=ai;
                        EVENTS(ai,2)=99;
                        EVENTS(ai,3)=0;
                        EVENTS(ai,4)=99;
                    end
                end
            end
        end
        
        % Make into BIDS format
        cd(rootdir);
        Bidsify_sub(subno,FileFormat,TASKNAME,EVENTS,EVENTSNames);
        
        eegdir=[rootdir,'sub-',subno,'\eeg\'];
        savedir=[rootdir,'derivatives\eegpreprocess\sub-',subno,'\eeg\'];
        saveCLEAN=[rootdir,'derivatives\processed\sub-',subno,'\eeg\'];
        behdir=[rootdir,'sub-',subno,'\beh\'];
        
        nbchannels = size(EEG.data,1);
        
        % Remove VEOG etc
        EEG.VEOG=squeeze(EEG.data(64,:,:));
        EEG.SCR=squeeze(EEG.data(65,:,:));
        EEG.EKG=squeeze(EEG.data(66,:,:));
        EEG.data=EEG.data(1:63,:,:);
        EEG.nbchan=63;
        EEG.chanlocs(66)=[]; EEG.chanlocs(65)=[]; EEG.chanlocs(64)=[];
        % Fix BV-specific issue - - - only needed for APPLE
        for ai=1:size(EEG.urevent,2), EEG.urevent(ai).bvtime=EEG.urevent(ai).bvmknum; end
        for ai=1:size(EEG.event,2), EEG.event(ai).bvtime=EEG.event(ai).bvmknum; end
        for ai=1:size(EEG.epoch,2), EEG.epoch(ai).eventbvtime=EEG.epoch(ai).eventbvmknum; end
        % Add CPz
        EEG = pop_chanedit(EEG,'append',63,'changefield',{64 'labels' 'CPz'});
        EEG = pop_chanedit(EEG,'lookup', locpath);
        % Re-Ref to Average Ref and recover CPz
        EEG = pop_reref(EEG,[],'refloc',struct('labels',{'CPz'},'type',{''},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},...
            'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{64},'ref',{''}),'keepref','on');
        % Remove everything else NOW that CPz has been reconstructed from the total
        EEG.MASTOIDS = squeeze(mean(EEG.data([10,21],:,:),1));
        EEG.data = EEG.data([1:4,6:9,11:20,22:26,28:64],:,:);
        EEG.nbchan=60;
        EEG.chanlocs(27)=[];  EEG.chanlocs(21)=[];   EEG.chanlocs(10)=[];   EEG.chanlocs(5)=[];  % Have to be in this order!
        % Should probably re-ref to average again now that the contaminated channels are gone
        EEG = pop_reref(EEG,[]);
        % Remove mean
        EEG = pop_rmbase(EEG,[],[]);
        
        % ----------------------
        % Setup APPLE to interp chans, reject epochs, & ID bad ICs.  Output will be Avg ref'd and ICA'd.
        session=1;
        eeg_chans=1:60;
        Do_ICA=1;
        subno=subjs(subji);
        ref_chan=36;  % Re-Ref to FCz   [WEIRD STEP, BUT THIS IS FOR FASTER, which is a part of APPLE]
        EEG = pop_reref(EEG,ref_chan,'keepref','on');
        BAAAD_Epochs = EVENTS(EVENTS(:,2)==99,1);
        cd(savedir);
        
        % Run APPLE (will re-ref data to avg ref)
        [EEG,EEG.bad_chans,EEG.bad_epochs,EEG.bad_ICAs,EEG.PV]=APPLE_SoS(EEG,eeg_chans,ref_chan,Do_ICA,subno,EEG.VEOG,session,BAAAD_Epochs);
        
        % Delete Bad Channels from EVENTS
        Baddie = find(double(cell2mat(EEG.bad_epochs(1,1)) == 1));
        Baddie = [Baddie;EVENTS(EVENTS(:,2)==99,1)];
        GOODS=ones(1,size(EVENTS,1));  GOODS(Baddie)=0;
        EVENTS_CLEAN=EVENTS(logical(GOODS),:);
        
        % Save
        save([savedir,num2str(subno),'_',TASKNAME,'_ready.mat'],'EEG','EVENTS', 'EVENTS_CLEAN');
        % ----------------------

    else
        eegdir=[rootdir,'sub-',num2str(subno),'\eeg\'];
        behdir=[rootdir,'sub-',num2str(subno),'\beh\'];
        savedir=[rootdir,'derivatives\eegpreprocess\sub-',num2str(subno),'\eeg\'];
        saveCLEAN=[rootdir,'derivatives\processed\sub-',num2str(subno),'\eeg\'];
        behdir=[rootdir,'sub-',num2str(subno),'\beh\'];
        
    end  % end exist clean.mat
    
    
    if  ~exist([saveCLEAN, num2str(subno),'_SoS_L.mat'])
        
        load([savedir,num2str(subno),'_',TASKNAME,'_ready.mat'],'EEG','EVENTS_CLEAN');
        eeglab redraw;
        EVENTS = EVENTS_CLEAN;
        
        if ischar(subno)
            subno = str2num(subno);
        end
        
        % Remove bad ICA
        
%         bad_ICAs_To_Remove=1;
        [NUM,TXT,RAW] = xlsread('SoSBadICA.xlsx');
        BadICA = NUM(find(subno==NUM(:,1)),2:end);
        for ICAi = size(BadICA,2):-1:1
            if ~isnan(BadICA(ICAi))
                bad_ICAs_To_Remove(ICAi)=BadICA(ICAi);
            end
        end
        clear ICAi BadICA
        EEG = pop_subcomp( EEG, bad_ICAs_To_Remove, 0);
        
        % Set up times for each
        tx=-2000:1000/EEG.srate:5998;
        TFB1=find(tx==-300);  TFB2=find(tx==-200);
        B1=find(tx==-200);  B2=find(tx==0);
        T1=find(tx==-500);  T2=find(tx==1500);
        tx2disp=-500:1000/EEG.srate:1500;
        
        
        for REFi=1:2
            if REFi==1, Label_Ref='V';
            elseif REFi==2, Label_Ref='L';
                X = [EEG.chanlocs.X]; Y = [EEG.chanlocs.Y]; Z = [EEG.chanlocs.Z];
                [EEG.data,G,H] = laplacian_perrinX(EEG.data,X,Y,Z,[],1e-6);
            end
            
            dims = size(EEG.data);
            % Shift the data from stim-locked to FB-locked
            for ai=1:dims(3)
                if EVENTS(ai,2) == 2 || EVENTS(ai,2) == 4
                    RTsamples(ai)=round((EVENTS(ai,3)+90)/2);  % /2 to make samples at 500 Hz
                    EEG.FB(:,:,ai)=[squeeze(EEG.data(:,RTsamples(ai):end,ai)),zeros(dims(1),RTsamples(ai)-1)];
                else
                    RTsamples(ai)=round((EVENTS(ai,3))/2);  % /2 to make samples at 500 Hz
                    EEG.FB(:,:,ai)=[squeeze(EEG.data(:,RTsamples(ai):end,ai)),zeros(dims(1),RTsamples(ai)-1)];
                end
            end
            
            % ------------------------------- TF
            
            % Setup data params
            for CTFi = 1:size(ElectrodesForTF,2)
                Channels(CTFi) = find(strcmpi(ElectrodesForTF{CTFi},{EEG.chanlocs.labels}));
            end
            
            for ChanNum = 1:size(Channels,2)
                
                chani=Channels(1,ChanNum);
                
                % Setup Wavelet Params
                num_freqs=50;
                frex=logspace(.01,1.7,num_freqs);
                s=logspace(log10(3),log10(10),num_freqs)./(2*pi*frex);
                t=-2:1/EEG.srate:2;
                
                % Define Convolution Parameters
                dims = size(EEG.data);
                n_wavelet = length(t);
                n_data = dims(2)*dims(3);
                n_convolution = n_wavelet+n_data-1;
                n_conv_pow2 = pow2(nextpow2(n_convolution));
                half_of_wavelet_size = (n_wavelet-1)/2;
                
                % get FFT of data
                EEG_CUE_fft = fft(reshape(EEG.data(chani,:,:),1,n_data),n_conv_pow2);
                EEG_FB_fft = fft(reshape(EEG.FB(chani,:,:),1,n_data),n_conv_pow2);
                
                for fi=1:num_freqs
                    
                    wavelet = fft( exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*(s(fi)^2))) , n_conv_pow2 );
                    
                    % convolution
                    EEG_conv = ifft(wavelet.*EEG_CUE_fft);
                    EEG_conv = EEG_conv(1:n_convolution);
                    EEG_conv = EEG_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    EEG_CUE_conv = reshape(EEG_conv,dims(2),dims(3));  clear EEG_conv;
                    
                    EEG_conv = ifft(wavelet.*EEG_FB_fft);
                    EEG_conv = EEG_conv(1:n_convolution);
                    EEG_conv = EEG_conv(half_of_wavelet_size+1:end-half_of_wavelet_size);
                    EEG_FB_conv = reshape(EEG_conv,dims(2),dims(3));  clear EEG_conv;
                    
                    % Get baseline as PRE-CUE!!
                    BASE = mean(mean( abs( EEG_CUE_conv(TFB1:TFB2,:) ) .^2,1),2);
                    
                    % Get power
                    temp_POWER_REW1 = mean(abs(  EEG_FB_conv(T1:T2,EVENTS(:,2)==1) ).^2 ,2);
                    temp_POWER_REW2 = mean(abs(  EEG_FB_conv(T1:T2,EVENTS(:,2)==2) ).^2 ,2);
                    temp_POWER_NULL1 = mean(abs(  EEG_FB_conv(T1:T2,EVENTS(:,2)==3) ).^2 ,2);
                    temp_POWER_NULL2 = mean(abs(  EEG_FB_conv(T1:T2,EVENTS(:,2)==4) ).^2 ,2);
                    
                    % dB correct power by base
                    POWER.FB(ChanNum,fi,:,1) = 10*(log10(temp_POWER_NULL2) - log10(repmat(BASE,size(temp_POWER_NULL2))));
                    POWER.FB(ChanNum,fi,:,2) = 10*(log10(temp_POWER_REW1) - log10(repmat(BASE,size(temp_POWER_REW1))));
                    POWER.FB(ChanNum,fi,:,3) = 10*(log10(temp_POWER_REW2) - log10(repmat(BASE,size(temp_POWER_REW2))));
                    POWER.FB(ChanNum,fi,:,4) = 10*(log10(temp_POWER_NULL1) - log10(repmat(BASE,size(temp_POWER_NULL1))));
                    
                    % Get ITPC by condi (different time windows)
                    ITPC.FB(ChanNum,fi,:,1) = abs(mean(exp(1i*(  angle(EEG_FB_conv(T1:T2,EVENTS(:,2)==4))  )),2));
                    ITPC.FB(ChanNum,fi,:,2) = abs(mean(exp(1i*(  angle(EEG_FB_conv(T1:T2,EVENTS(:,2)==1))  )),2));
                    ITPC.FB(ChanNum,fi,:,3) = abs(mean(exp(1i*(  angle(EEG_FB_conv(T1:T2,EVENTS(:,2)==2))  )),2));
                    ITPC.FB(ChanNum,fi,:,4) = abs(mean(exp(1i*(  angle(EEG_FB_conv(T1:T2,EVENTS(:,2)==3))  )),2));
                    
                    % ----------------
                    clear temp* EEG_FB_conv wavelet BASE;
                    
                end
            end
            cd(saveCLEAN);
            figure; 
            subplot(2,1,1); hold on
            imagesc(tx2disp,[], squeeze(POWER.FB(1,:,:,2)-POWER.FB(1,:,:,1)) );
            plot([0 0],[1 50],'k:')
            set(gca,'clim',[-4 4],'xlim',[-500,1000],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
            title(['Sight REW - NULL    Sx:',num2str(subno)]); cbar
            subplot(2,1,2); hold on
            imagesc(tx2disp,[], squeeze(POWER.FB(1,:,:,4)-POWER.FB(1,:,:,3)) );
            plot([0 0],[1 50],'k:')
            set(gca,'clim',[-4 4],'xlim',[-500,1000],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
            title(['Sound REW - NULL    Sx:',num2str(subno)]); cbar
            
            saveas(gcf, [num2str(subno),'_SnS_FB_TF',Label_Ref,'.png'],'png');
            
            % ------------------------------- ERPs
            
            % Let's just do this for display
            dims=size(EEG.data);
            EEG.FILT=eegfilt(EEG.FB,500,[],20);
            EEG.FILT=reshape(EEG.FILT,dims(1),dims(2),dims(3));
            
            % Basecor your ERPs here so they are pretty.
            BASE_cue=squeeze(  mean(EEG.FILT(:,B1:B2,:),2)  );
            
            for ai=1:dims(1)
                EEG.FILT(ai,:,:)=squeeze(EEG.FILT(ai,:,:))-repmat( BASE_cue(ai,:),dims(2),1 );
            end
            
            % AVERAGE
            ERPs.FBbyCUE(:,:,1)=squeeze(mean(  EEG.FILT(:,T1:T2,EVENTS(:,2)==4) ,3));
            ERPs.FBbyCUE(:,:,2)=squeeze(mean(  EEG.FILT(:,T1:T2,EVENTS(:,2)==1) ,3));
            ERPs.FBbyCUE(:,:,3)=squeeze(mean(  EEG.FILT(:,T1:T2,EVENTS(:,2)==2) ,3));
            ERPs.FBbyCUE(:,:,4)=squeeze(mean(  EEG.FILT(:,T1:T2,EVENTS(:,2)==3) ,3));
            
            % SUPER IMPORTANT KINDA HOUSEKEEPING AT LEAST VERY GOOD PRACTICE
            EEG.FILT=[];
            
            
            % -----------------------------------
            chani=36;
            FBCOLS = {'r','g'};
            figure; 
            subplot(2,1,1); hold on
            for condi=1:2
                plot(tx2disp,squeeze(ERPs.FBbyCUE(chani,:,condi)),FBCOLS{condi},'linewidth',2);
            end
            title('Sight FB by cue');
            subplot(2,1,2); hold on
            for condi=1:2
                plot(tx2disp,squeeze(ERPs.FBbyCUE(chani,:,condi+2)),[FBCOLS{condi},':'],'linewidth',2);
            end
            title('Sight FB by cue');
            
            saveas(gcf, [num2str(subno),'_SoS_FB_ERP_',Label_Ref,'.png'],'png');
            
            close all;
            
            
            save([saveCLEAN, num2str(subno),'_SoS_',Label_Ref,'.mat'],'EVENTS','ERPs','POWER','ITPC');
            clearvars -except ElectrodesForTF TaskLabel TASKNAME FileFormat DATE rootdir USEsubjs sourcedir MEGA_BEH ALL CB1 CB2 homedir subdir locpath eegdir raweegdir savedir saveCLEAN sx si subjs subji subno EEG EVENTS tx TFB1 TFB2 B1 B2 T1 T2 tx2disp FILTERS USE;
        end % End Voltage vs. Laplacian
        
        
    end  % end TF/ERP
    
    
    load([saveCLEAN, num2str(subno),'_SoS_V.mat'],'ERPs','POWER','ITPC');
    
    if ischar(subno)
        subno = str2num(subno);
    end
    
    USEsubjs(ALL) = subjs(subji);
    
    MEGA_ERPs.FBbyCUE(ALL,:,:,:) = ERPs.FBbyCUE;
    MEGA_TF.FB(ALL,:,:,:,:) = POWER.FB;
    MEGA_ITPC.FB(ALL,:,:,:,:) = ITPC.FB;
    clear ERPS POWER ITPC
    
    load([saveCLEAN, num2str(subno),'_SoS_L.mat'],'ERPs','POWER','ITPC');
    
    MEGA_ERPs.LAP_FBbyCUE(ALL,:,:,:) = ERPs.FBbyCUE;
    MEGA_TF.LAP_FB(ALL,:,:,:,:) = POWER.FB;
    MEGA_ITPC.LAP_FB(ALL,:,:,:,:) = ITPC.FB;
    ALL = ALL+1;
    
    
    clearvars -except ElectrodesForTF TaskLabel TASKNAME FileFormat DATE rootdir USEsubjs MEGA_BEH sourcedir homedir subdir locpath sx si subjs subji subno ...
        MEGA_ERPs MEGA_TF MEGA_ITPC UP_ERPs UP_TF UP_ITPC MED_ERPs MED_TF ...
        MED_ITPC CB1 CB2 ALL USE
end

clear USE


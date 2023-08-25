function [Model_AIC] = TestRewPModel(params,model,Pun_ERP,Rew_ERP,tx2disp,srate,n,time1,time2,UseBaseline); 

if strcmp(model,'Full')
    %% ######## Create full cosine wave
    
    % Pull Params
    FREQ_fit=params(1);
    SCALE_fit=params(2);
    TSHIFT_fit=round(params(3));
    NumOfParams = 3;
    NumOfSamples = (time2-time1)./2;
    
    % Get time
    full_ti= -((1/FREQ_fit)*.25):(1000/srate)/1000:((1/FREQ_fit)*.75);

    % Fit Freq
    full_burst = cos(2*pi*FREQ_fit.*full_ti);

    % Amplitude scaling
    SCALE_Full = full_burst.*SCALE_fit;

    % Add Baseline and remove tail if beyond time window
    AddBaseline  = find(tx2disp(1,:)==UseBaseline)+TSHIFT_fit;
    TIME_Full = [zeros(1,AddBaseline),SCALE_Full,zeros(1,(751-(length(SCALE_Full)+AddBaseline)))];
    TIME_Full = TIME_Full(:,1:751);
    
    % Add to Punishment condition
    SimRewP_Full = Pun_ERP(n,:) + TIME_Full;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_Full = sum((Rew_ERP(n,T1:T2)-SimRewP_Full(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_Full./NumOfSamples) + 2.* NumOfParams;


elseif strcmp(model,'Half')
    %% ######## Create half cosine wave
    
    % Pull Params
    FREQ_fit=params(1);
    SCALE_fit=params(2);
    TSHIFT_fit=round(params(3));
    NumOfParams = 3;
    NumOfSamples = (time2-time1)./2;
    
    % Get time
    half_ti= -((1/FREQ_fit)*.25):(1000/srate)/1000:((1/FREQ_fit)*.25);
        
    % Fit Freq
    half_burst = cos(2*pi*FREQ_fit.*half_ti);
    
    % Amplitude scaling
    SCALE_Half = half_burst.*SCALE_fit;
    
    % Add Baseline and remove tail if beyond time window
    AddBaseline  = find(tx2disp(1,:)==UseBaseline)+TSHIFT_fit;
    TIME_Half = [zeros(1,AddBaseline),SCALE_Half,zeros(1,(751-(length(SCALE_Half)+AddBaseline)))];
    TIME_Half = TIME_Half(:,1:751);
    
    % Add to Punishment condition
    SimRewP_Half = Pun_ERP(n,:) + TIME_Half;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_Half = sum((Rew_ERP(n,T1:T2)-SimRewP_Half(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_Half./NumOfSamples) + 2.* NumOfParams;
        
elseif strcmp(model,'Gauss')
    %% ######## Create Gaussian
    
    % Pull Params
    FREQ_fit=params(1);
    SCALE_fit=params(2);
    TSHIFT_fit=round(params(3));
    NumOfParams = 3;
    NumOfSamples = (time2-time1)./2;
    Sigma = 2;
    
    % Get time
    Gauss_ti = -5:.001:5;
        
    % Fit Freq
    s=Sigma./(2*pi.*FREQ_fit);
    Gauss=exp(-Gauss_ti.^2./(2*s^2));
    Gauss=Gauss(1:4:end); % downsample
    
    % Amplitude scaling
    SCALE_Gauss=Gauss.*SCALE_fit;
    if SCALE_fit >=0
        SCALE_Gauss=SCALE_Gauss(1,SCALE_Gauss(1,:)>0.001); % round tails for movement
    elseif SCALE_fit <0
        SCALE_Gauss=SCALE_Gauss(1,SCALE_Gauss(1,:)<(-0.001));
    end
    
    
    % Add Baseline and remove tail if beyond time window
    AddBaseline  = TSHIFT_fit;
    TIME_Gauss = [zeros(1,AddBaseline),SCALE_Gauss,zeros(1,(751-(length(SCALE_Gauss)+AddBaseline)))];
    TIME_Gauss = TIME_Gauss(:,1:751);
    
    % Add to Punishment condition
    SimRewP_Gauss = Pun_ERP(n,:) + TIME_Gauss;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_Gauss = sum((Rew_ERP(n,T1:T2)-SimRewP_Gauss(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_Gauss./NumOfSamples) + 2.* NumOfParams;
    
elseif strcmp(model,'Wavelet')

    %% ######## Create Wavelet
    
    % Pull Params
    FREQ_fit=params(1);
    SCALE_fit=params(2);
    TAPER_fit=round(params(3));
    TSHIFT_fit=round(params(4));
    NumOfParams = 4;
    NumOfSamples = (time2-time1)./2;
    
    % Get time
    Wavelet_ti=-1:.001:1;
    
    % Fit Freq
    Delta=exp(2*1i*pi*FREQ_fit.*Wavelet_ti);
    s=2./(2*pi.*FREQ_fit);
    Delta_Taper=exp(-Wavelet_ti.^2./(2*s^2));
    
    % Fit Taper over the sine
    Delta_Taper_Shift=[zeros(1,TAPER_fit),Delta_Taper(1:end-TAPER_fit)];
    Deltalet=Delta.*Delta_Taper_Shift;
    
    % Manage temporal duration
    Deltalet=real(Deltalet(500:1500));    % 1001 samples
    % Re-sample based on srate of data
    Deltalet=resample(Deltalet,1,1000/srate);    
    
    % Amplitude scaling
    Deltalet_scale=Deltalet.*SCALE_fit;
    Deltalet_scale=round(Deltalet_scale,3); % cut off tails to make fit, technically imprecise
    Deltalet_scale(:,Deltalet_scale(:,:)==0)=[]; % cut off tails to make fit, technically imprecise
    
    % Add Baseline and remove tail if beyond time window
    AddBaseline  = TSHIFT_fit;
    TIME_Wavelet = [zeros(1,AddBaseline),Deltalet_scale,zeros(1,(751-(length(Deltalet_scale)+AddBaseline)))];
    TIME_Wavelet = TIME_Wavelet(:,1:751);
    
    % Add to Punishment condition
    SimRewP_Wavelet = Pun_ERP(n,:) + TIME_Wavelet;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_Wavelet = sum((Rew_ERP(n,T1:T2)-SimRewP_Wavelet(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_Wavelet./NumOfSamples) + 2.* NumOfParams;

elseif strcmp(model,'Gamma')
    %% ######## Create Gamma
    
    % Pull Params
    SHAPE_fit=params(1);
    RATE_fit=params(2);
    SCALE_fit=params(3);
    TSHIFT_fit=round(params(4));
    NumOfParams = 4;
    NumOfSamples = (time2-time1)./2;
    
    % Get time
    Gamma_ti = 0:.001:100;
        
    % Fit shape and rate params
    Gamma = makedist('Gamma','a',SHAPE_fit,'b',RATE_fit);
    Gamma_Dist=pdf(Gamma,Gamma_ti);
    Gamma_Dist=Gamma_Dist(1:60:end); % downsample to fit with ERPs
    Gamma_Dist=Gamma_Dist(1,Gamma_Dist(1,:)>0.001);
    
    % Amplitude scaling
    SCALE_Gamma=Gamma_Dist.*SCALE_fit;
    
    % Add Baseline and remove tail if beyond time window
    AddBaseline  = find(tx2disp(1,:)==UseBaseline)+TSHIFT_fit;
    TIME_Gamma = [zeros(1,AddBaseline),SCALE_Gamma,zeros(1,(751-(length(SCALE_Gamma)+AddBaseline)))];
    TIME_Gamma = TIME_Gamma(:,1:751);
    
    % Add to Punishment condition
    SimRewP_Gamma = Pun_ERP(n,:) + TIME_Gamma;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_Gamma = sum((Rew_ERP(n,T1:T2)-SimRewP_Gamma(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_Gamma./NumOfSamples) + 2.* NumOfParams;
    
elseif strcmp(model,'TwoHalf')
    %% ######## Create half cosine wave
    
    % Pull Params
    FREQ_fit=params(1);
    SCALE_fit=params(2);
    TSHIFT_fit=round(params(3));
    FREQ2_fit=params(4);
    SCALE2_fit=params(5);
    TSHIFT2_fit=round(params(6));
    NumOfParams = 6;
    NumOfSamples = (time2-time1)./2;
    
    % Get time
    Firsthalf_ti= -((1/FREQ_fit)*.25):(1000/srate)/1000:((1/FREQ_fit)*.25);
        
    % Fit Freq
    Firsthalf_burst = cos(2*pi*FREQ_fit.*Firsthalf_ti);
    
    % Amplitude scaling
    FirstSCALE_Half = Firsthalf_burst.*SCALE_fit;
    
    % Add Baseline and remove tail if beyond time window
    FirstAddBaseline  = find(tx2disp(1,:)==UseBaseline)+TSHIFT_fit;
    FirstTIME_Half = [zeros(1,FirstAddBaseline),FirstSCALE_Half,zeros(1,(751-(length(FirstSCALE_Half)+FirstAddBaseline)))];
    FirstTIME_Half = FirstTIME_Half(:,1:751);
    
    % Get Second time
    Secondhalf_ti= -((1/FREQ2_fit)*.25):(1000/srate)/1000:((1/FREQ2_fit)*.25);
        
    % Fit Freq
    Secondhalf_burst = cos(2*pi*FREQ2_fit.*Secondhalf_ti);
    
    % Amplitude scaling
    SecondSCALE_Half = -1.*Secondhalf_burst.*SCALE2_fit;
    
    % Add Baseline and remove tail if beyond time window
    SecondAddBaseline  = FirstAddBaseline+TSHIFT2_fit;
    SecondTIME_Half = [zeros(1,SecondAddBaseline),SecondSCALE_Half,zeros(1,(751-(length(SecondSCALE_Half)+SecondAddBaseline)))];
    SecondTIME_Half = SecondTIME_Half(:,1:751);
    
    TIME_2Half = FirstTIME_Half + SecondTIME_Half;
    
    % Add to Punishment condition
    SimRewP_2Half = Pun_ERP(n,:) + TIME_2Half;
    
    % Compare simulated with Real
    T1 = find(tx2disp(1,:)==time1); T2 = find(tx2disp(1,:)==time2);
    SSE_2Half = sum((Rew_ERP(n,T1:T2)-SimRewP_2Half(1,T1:T2)).^2);
    Model_AIC = NumOfSamples .* log(SSE_2Half./NumOfSamples) + 2.* NumOfParams;
        
    
end


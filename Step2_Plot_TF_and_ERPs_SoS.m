%% ERPs 
tx2disp=-500:2:1000;

Channels = [2 36 21];   % [2 36 21 60 11 45 15]; % Fz FCz Cz CPz Pz POz Oz
T1 = 100; T2 = 250;
figure;

clear TitleTextFB
TitleTextFB = 'Reward vs Null';
LineColor = ['b'; 'k'; 'g'; 'r'];
Red = [1 0.75 0.8]; Blue = [0.529 0.808 0.98]; Green = [0.59 0.98 0.59]; Black = [0.862 0.862 0.862];
SDLineColor = [Blue; Black; Green; Red];
BigN = size(MEGA_ERPs.FBbyCUE,1);

    
subplot(4,1,1, 'Position',[0.5 0.9 0.01 0.01]);
title('Reward (Sight = blue, Sound = green) vs. Null (Sight = black, Sound = red)');
axis off;

for chani = 1:3
    
    if chani == 1
        clear TitleTextCh
        TitleTextCh = 'Fz';
        PlotSpace = [4;7;10];
    elseif chani == 2
        clear TitleTextCh
        TitleTextCh = 'FCz';
        PlotSpace = [5;8;11];
    elseif chani == 3
        clear TitleTextCh
        TitleTextCh = 'Cz';
        PlotSpace = [6;9;12];
        
    end
    
    SiRew = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,2);
    SiPun = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,4);
    SoRew = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,3);
    SoPun = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,1);
    
    Time1SD = squeeze(  std( SiRew  ,0,1) ./ sqrt(BigN)  )';
    Time2SD = squeeze(  std( SiPun  ,0,1) ./ sqrt(BigN)  )';
    Time3SD = squeeze(  std( SoRew  ,0,1) ./ sqrt(BigN)  )';
    Time4SD = squeeze(  std( SoPun  ,0,1) ./ sqrt(BigN)  )';
    
    subplot(4,1,chani+1); hold on
    
    rectangle('Position',[T1,-10,(T2-T1),20],'Curvature',0.1,'FaceColor',[.5 .5 .5])
    
    xFill = [tx2disp fliplr(tx2disp)];
    yFill = [(squeeze(   mean(  SiRew  ,1))'+Time1SD) fliplr((squeeze(   mean(  SiRew  ,1))'-Time1SD))];
    h = fill(xFill,yFill,SDLineColor(1,:),'edgecolor','none');
    set(h,'facealpha',0.5)
    yFill = [(squeeze(   mean(  SiPun  ,1))'+Time2SD) fliplr((squeeze(   mean(  SiPun  ,1))'-Time2SD))];
    h = fill(xFill,yFill,SDLineColor(2,:),'edgecolor','none');
    set(h,'facealpha',0.5)
    xFill = [tx2disp fliplr(tx2disp)];
    yFill = [(squeeze(   mean(  SoRew  ,1))'+Time3SD) fliplr((squeeze(   mean(  SoRew  ,1))'-Time3SD))];
    h = fill(xFill,yFill,SDLineColor(3,:),'edgecolor','none');
    set(h,'facealpha',0.5)
    yFill = [(squeeze(   mean(  SoPun  ,1))'+Time4SD) fliplr((squeeze(   mean(  SoPun  ,1))'-Time4SD))];
    h = fill(xFill,yFill,SDLineColor(4,:),'edgecolor','none');
    set(h,'facealpha',0.5)
    
    plot(-500:1000/500:1000,squeeze(  mean( SiRew ,1)  ),LineColor(1,:),'linewidth',2);
    plot(-500:1000/500:1000,squeeze(  mean( SiPun ,1)  ),LineColor(2,:),'linewidth',2);
    plot(-500:1000/500:1000,squeeze(  mean( SoRew ,1)  ),LineColor(3,:),'linewidth',2);
    plot(-500:1000/500:1000,squeeze(  mean( SoPun ,1)  ),LineColor(4,:),'linewidth',2);
    plot([0 0],[-6 6],'k:');
    plot([-500 1000], [0 0], 'k:');
    plot([-200 -200],[-6 6], 'b:');
    set(gca,'xlim',[-500 1000],'ylim',[-10 10],'xtick',[-500:500:1500])
    pcrit=.05;
%     [h,P,ci,STATS]=ttest( SiRew,SiPun );
%     P(P>pcrit)=NaN;  P(P<=pcrit)=1;
%     plot(tx2disp,(-7 .* squeeze(P)),'k','linewidth',3);  clear H P CI STATS
%     [h,P,ci,STATS]=ttest( SoRew,SoPun );
%     P(P>pcrit)=NaN;  P(P<=pcrit)=1;
%     plot(tx2disp,(-9 .* squeeze(P)),'r','linewidth',3);  clear H P CI STATS
    title([TitleTextCh]);
    clear Time1SD Time2SD Time3SD Time4SD ONE TWO THREE FOUR
    
end

% [1:4 7:10 13:16 19:22]
% cd(['C:\Users\GA217C\Desktop\Projects\GeneProject\Pics\']);
% saveas(gcf,['NSF_ERP_',TitleTextFB,'_CoMT'],'png')
    

%% TF
figure; hold on
frex=logspace(.01,1.7,50);
ttesttype='within';

SiRew = MEGA_TF.FB(:,2,:,:,2);
SiPun = MEGA_TF.FB(:,2,:,:,4);
SoRew = MEGA_TF.FB(:,2,:,:,3);
SoPun = MEGA_TF.FB(:,2,:,:,1);

subplot(2,3,1); hold on
imagesc(-500:2:1500,[],squeeze(   mean(SiRew,1)   )); axis xy
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title('Sight Rew');

subplot(2,3,2); hold on
imagesc(-500:2:1500,[],squeeze(   mean(SiPun,1)   )); axis xy
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title('Sight Null');

subplot(2,3,4); hold on
imagesc(-500:2:1500,[],squeeze(   mean(SoRew,1)   )); axis xy
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title('Sound Rew');

subplot(2,3,5); hold on
imagesc(-500:2:1500,[],squeeze(   mean(SoPun,1)   )); axis xy
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title('Sound Null');
        
subplot(2,3,3); hold on
imagesc(-500:2:1500,[],squeeze( mean(SiRew,1) ) - squeeze( mean(SiPun,1) )); axis xy
[H,P,CI,STATS]=ttest2(squeeze(SiRew),squeeze(SiPun));
A=squeeze(SiRew);
B=squeeze(SiPun);
[Corrected_P] = Run_Thresh_2D(A,B,ttesttype); clear A B
contour(-500:2:1500,1:50,Corrected_P,'k','linewidth',2); clear Corrected_P H P CI STATS;
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title(['Sight Rew - Null']); cbar

subplot(2,3,6); hold on
imagesc(-500:2:1500,[],squeeze( mean(SoRew,1) ) - squeeze( mean(SoPun,1) )); axis xy
[H,P,CI,STATS]=ttest2(squeeze(SoRew),squeeze(SoPun));
A=squeeze(SoRew);
B=squeeze(SoPun);
[Corrected_P] = Run_Thresh_2D(A,B,ttesttype); clear A B
contour(-500:2:1500,1:50,Corrected_P,'k','linewidth',2); clear Corrected_P H P CI STATS;
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1500],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title(['Sound Rew - Null']); cbar



%% Topo's
tx2disp=-500:2:1500;
load BV_ChanLocs_60

% Sight ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SightT1 = 250; SightT2 = 400;
Sitopo1=find(tx2disp==SightT1); Sitopo2=find(tx2disp==SightT2);
% Sound ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SoundT1 = 100; SoundT2 = 250;
Sotopo1=find(tx2disp==SoundT1); Sotopo2=find(tx2disp==SoundT2);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get ERP's
SiRew = MEGA_ERPs.FBbyCUE(:,:,Sitopo1:Sitopo2,2);
SiPun = MEGA_ERPs.FBbyCUE(:,:,Sitopo1:Sitopo2,4);
SoRew = MEGA_ERPs.FBbyCUE(:,:,Sotopo1:Sotopo2,3);
SoPun = MEGA_ERPs.FBbyCUE(:,:,Sotopo1:Sotopo2,1);

SET1=squeeze(  mean(SiRew,3)  ); SET2=squeeze(  mean(SiPun,3)  );
SET3=squeeze(  mean(SoRew,3)  ); SET4=squeeze(  mean(SoPun,3)  );
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot

figure;
subplot(2,3,1);
topoplot(mean(SET1,1),BV_Chanlocs_60,'maplimits',[-3 3]);

subplot(2,3,2);
topoplot(mean(SET2,1),BV_Chanlocs_60,'maplimits',[-3 3]);

title(['Sight Rew - Null']);
DATA=squeeze(mean(SET1,1)-mean(SET2,1));  P=ones(1,60);
[H,P,CI,STATS]=ttest(SET1,SET2); P(P>=.05)=NaN; P(P<.05)=1; P(isnan(P))=0;
subplot(2,3,3);
topoplot(DATA,BV_Chanlocs_60,'maplimits',[-3 3],'emarker2', {find(P),'d','k',10,1});
cbar;

subplot(2,3,4);
topoplot(mean(SET3,1),BV_Chanlocs_60,'maplimits',[-3 3]);

subplot(2,3,5);
topoplot(mean(SET4,1),BV_Chanlocs_60,'maplimits',[-3 3]);

title(['Sound Rew - Null']);
DATA=squeeze(mean(SET3,1)-mean(SET4,1));  P=ones(1,60);
[H,P,CI,STATS]=ttest(SET3,SET4); P(P>=.05)=NaN; P(P<.05)=1; P(isnan(P))=0;
subplot(2,3,6);
topoplot(DATA,BV_Chanlocs_60,'maplimits',[-3 3],'emarker2', {find(P),'d','k',10,1});
cbar;

%% Combined Figures
tx2disp=-500:2:1000;

chani = 3;

if chani == 1
    ChannelText = 'Fz';
elseif chani == 2
    ChannelText = 'FCz';
elseif chani == 3
    ChannelText = 'Cz';
end

% Sight Stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SiRew_ERP = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,2);
SiPun_ERP = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,4);
Time1SD = squeeze(  std( SiRew_ERP  ,0,1) ./ sqrt(BigN)  )';
Time2SD = squeeze(  std( SiPun_ERP  ,0,1) ./ sqrt(BigN)  )';
SiRew_TOPO = MEGA_ERPs.FBbyCUE(:,:,Sitopo1:Sitopo2,2);
SiPun_TOPO = MEGA_ERPs.FBbyCUE(:,:,Sitopo1:Sitopo2,4);
SiRew_TF = MEGA_TF.FB(:,3,:,1:751,2);
SiPun_TF = MEGA_TF.FB(:,3,:,1:751,4);    
SET1=squeeze(  mean(SiRew_TOPO,3)  ); SET2=squeeze(  mean(SiPun_TOPO,3)  );

% Sound Stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SoRew_ERP = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,3);
SoPun_ERP = MEGA_ERPs.FBbyCUE(:,Channels(chani),1:751,1);
Time3SD = squeeze(  std( SoRew_ERP  ,0,1) ./ sqrt(BigN)  )';
Time4SD = squeeze(  std( SoPun_ERP  ,0,1) ./ sqrt(BigN)  )';
SoRew_TOPO = MEGA_ERPs.FBbyCUE(:,:,Sotopo1:Sotopo2,3);
SoPun_TOPO = MEGA_ERPs.FBbyCUE(:,:,Sotopo1:Sotopo2,1);
SoRew_TF = MEGA_TF.FB(:,3,:,1:751,3);
SoPun_TF = MEGA_TF.FB(:,3,:,1:751,1);
SET3=squeeze(  mean(SoRew_TOPO,3)  ); SET4=squeeze(  mean(SoPun_TOPO,3)  );

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure;
subplot(2,4,1:2); hold on

rectangle('Position',[T2,-10,(T2-T1),20],'Curvature',0.1,'FaceColor',[.5 .5 .5])

xFill = [tx2disp fliplr(tx2disp)];
yFill = [(squeeze(   mean(  SiRew_ERP  ,1))'+Time1SD) fliplr((squeeze(   mean(  SiRew_ERP  ,1))'-Time1SD))];
h = fill(xFill,yFill,SDLineColor(3,:),'edgecolor','none');
set(h,'facealpha',0.5)
yFill = [(squeeze(   mean(  SiPun_ERP  ,1))'+Time2SD) fliplr((squeeze(   mean(  SiPun_ERP  ,1))'-Time2SD))];
h = fill(xFill,yFill,SDLineColor(4,:),'edgecolor','none');
set(h,'facealpha',0.5)
plot(-500:1000/500:1000,squeeze(  mean( SiRew_ERP ,1)  ),LineColor(3,:),'linewidth',2);
plot(-500:1000/500:1000,squeeze(  mean( SiPun_ERP ,1)  ),LineColor(4,:),'linewidth',2);
plot(-500:1000/500:1000,squeeze(  mean( SiRew_ERP ,1) - mean( SiPun_ERP ,1)  ),'b','linewidth',2);
plot([0 0],[-6 6],'k:');
plot([-500 1000], [0 0], 'k:');
plot([-200 -200],[-6 6], 'b:');
set(gca,'xlim',[-500 1000],'ylim',[-10 10],'xtick',[-500:200:1000])
pcrit=.05;
[h,P,ci,STATS]=ttest( SiRew_ERP,SiPun_ERP );
P(P>pcrit)=NaN;  P(P<=pcrit)=1;
title(['Sight ',ChannelText]);
clear P

subplot(2,4,5:6); hold on

rectangle('Position',[T1,-10,(T2-T1),20],'Curvature',0.1,'FaceColor',[.5 .5 .5])

xFill = [tx2disp fliplr(tx2disp)];
yFill = [(squeeze(   mean(  SoRew_ERP  ,1))'+Time1SD) fliplr((squeeze(   mean(  SoRew_ERP  ,1))'-Time1SD))];
h = fill(xFill,yFill,SDLineColor(3,:),'edgecolor','none');
set(h,'facealpha',0.5)
yFill = [(squeeze(   mean(  SoPun_ERP  ,1))'+Time2SD) fliplr((squeeze(   mean(  SoPun_ERP  ,1))'-Time2SD))];
h = fill(xFill,yFill,SDLineColor(4,:),'edgecolor','none');
set(h,'facealpha',0.5)
plot(-500:1000/500:1000,squeeze(  mean( SoRew_ERP ,1)  ),LineColor(3,:),'linewidth',2);
plot(-500:1000/500:1000,squeeze(  mean( SoPun_ERP ,1)  ),LineColor(4,:),'linewidth',2);
plot(-500:1000/500:1000,squeeze(  mean( SoRew_ERP ,1) - mean( SoPun_ERP ,1)  ),'b','linewidth',2);
plot([0 0],[-6 6],'k:');
plot([-500 1000], [0 0], 'k:');
plot([-200 -200],[-6 6], 'b:');
set(gca,'xlim',[-500 1000],'ylim',[-10 10],'xtick',[-500:200:1000])
pcrit=.05;
[h,P,ci,STATS]=ttest( SoRew_ERP,SoPun_ERP );
P(P>pcrit)=NaN;  P(P<=pcrit)=1;
title(['Sound ',ChannelText]);
clear P

subplot(2,4,4); hold on
imagesc(-500:2:1000,[],squeeze( mean(SiRew_TF,1) ) - squeeze( mean(SiPun_TF,1) )); axis xy
[H,P,CI,STATS]=ttest2(squeeze(SiRew_TF),squeeze(SiPun_TF));
A=squeeze(SiRew_TF);
B=squeeze(SiPun_TF);
[Corrected_P] = Run_Thresh_2D_TREVOR(A,B,BigN,BigN,ttesttype); clear A B
contour(-500:2:1000,1:50,Corrected_P,'k','linewidth',2); clear Corrected_P H P CI STATS;
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1000],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title([ChannelText,' Sight Rew - Null']); cbar

subplot(2,4,8); hold on
imagesc(-500:2:1000,[],squeeze( mean(SoRew_TF,1) ) - squeeze( mean(SoPun_TF,1) )); axis xy
[H,P,CI,STATS]=ttest2(squeeze(SoRew_TF),squeeze(SoPun_TF));
A=squeeze(SoRew_TF);
B=squeeze(SoPun_TF);
[Corrected_P] = Run_Thresh_2D_TREVOR(A,B,BigN,BigN,ttesttype); clear A B
contour(-500:2:1000,1:50,Corrected_P,'k','linewidth',2); clear Corrected_P H P CI STATS;
plot([0 0],[1 50],'k:')
set(gca,'clim',[-3 3],'xlim',[-500,1000],'ylim',[1 50],'YTick',1:4:length(frex),'YTickLabel',round(frex(1:4:end)));
title([ChannelText,' Sound Rew - Null']); cbar


DATA=squeeze(mean(SET1,1)-mean(SET2,1));  P=ones(1,60);
[H,P,CI,STATS]=ttest(SET1,SET2); P(P>=.05)=NaN; P(P<.05)=1; P(isnan(P))=0;
subplot(2,4,3);
topoplot(DATA,BV_Chanlocs_60,'maplimits',[-3 3],'emarker2', {find(P),'d','k',10,1});
title(['Sight Rew - Null']); cbar;

DATA=squeeze(mean(SET3,1)-mean(SET4,1));  P=ones(1,60);
[H,P,CI,STATS]=ttest(SET3,SET4); P(P>=.05)=NaN; P(P<.05)=1; P(isnan(P))=0;
subplot(2,4,7);
topoplot(DATA,BV_Chanlocs_60,'maplimits',[-3 3],'emarker2', {find(P),'d','k',10,1});
title(['Sound Rew - Null']); cbar;


%%
tx2disp=-500:2:1000;

Channels = [2 36 21];   % [2 36 21 60 11 45 15]; % Fz FCz Cz CPz Pz POz Oz
T1 = 100; T2 = 250;

clear TitleTextFB
TitleTextFB = 'Reward vs Null';
LineColor = ['b'; 'k'; 'g'; 'r'];
Red = [1 0.75 0.8]; Blue = [0.529 0.808 0.98]; Green = [0.59 0.98 0.59]; Black = [0.862 0.862 0.862];
SDLineColor = [Blue; Black; Green; Red];
BigN = size(MEGA_ERPs.FBbyCUE,1);

    
subplot(4,1,1, 'Position',[0.5 0.9 0.01 0.01]);
title('Reward (Sight = blue, Sound = green) vs. Null (Sight = black, Sound = red)');
axis off;

for chani = 1:3
    
figure;

    if chani == 1
        clear TitleTextCh
        TitleTextCh = 'Fz';
        PlotSpace = [4;7;10];
    elseif chani == 2
        clear TitleTextCh
        TitleTextCh = 'FCz';
        PlotSpace = [5;8;11];
    elseif chani == 3
        clear TitleTextCh
        TitleTextCh = 'Cz';
        PlotSpace = [6;9;12];
        
    end
    
    title(TitleTextCh);

    
    for subjectnumbah = 1:22
    SiRew = MEGA_ERPs.FBbyCUE(subjectnumbah,Channels(chani),1:751,2);
    SiPun = MEGA_ERPs.FBbyCUE(subjectnumbah,Channels(chani),1:751,4);
    SoRew = MEGA_ERPs.FBbyCUE(subjectnumbah,Channels(chani),1:751,3);
    SoPun = MEGA_ERPs.FBbyCUE(subjectnumbah,Channels(chani),1:751,1);
    
    subplot(2,2,1); hold on
    plot(-500:1000/500:1000,squeeze(  mean( SiRew ,1)  ),LineColor(1,:),'linewidth',2);
    plot([0 0],[-6 6],'k:');
    plot([-500 1000], [0 0], 'k:');
    plot([-200 -200],[-6 6], 'b:');
    set(gca,'xlim',[-500 1000],'ylim',[-20 20],'xtick',[-500:500:1500])
    
    subplot(2,2,2); hold on
    plot(-500:1000/500:1000,squeeze(  mean( SiPun ,1)  ),LineColor(2,:),'linewidth',2);
    plot([0 0],[-6 6],'k:');
    plot([-500 1000], [0 0], 'k:');
    plot([-200 -200],[-6 6], 'b:');
    set(gca,'xlim',[-500 1000],'ylim',[-20 20],'xtick',[-500:500:1500])
    
    subplot(2,2,3); hold on
    plot(-500:1000/500:1000,squeeze(  mean( SoRew ,1)  ),LineColor(3,:),'linewidth',2);
    plot([0 0],[-6 6],'k:');
    plot([-500 1000], [0 0], 'k:');
    plot([-200 -200],[-6 6], 'b:');
    set(gca,'xlim',[-500 1000],'ylim',[-20 20],'xtick',[-500:500:1500])
    
    subplot(2,2,4); hold on
    plot(-500:1000/500:1000,squeeze(  mean( SoPun ,1)  ),LineColor(4,:),'linewidth',2);
    plot([0 0],[-6 6],'k:');
    plot([-500 1000], [0 0], 'k:');
    plot([-200 -200],[-6 6], 'b:');
    set(gca,'xlim',[-500 1000],'ylim',[-20 20],'xtick',[-500:500:1500])
    
    end
    
end


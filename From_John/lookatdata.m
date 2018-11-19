% Load data and visual rejection

%addpath ~/Dropbox/Matlab/library/misc/
%addpath ~/Dropbox/Matlab/thirdparty/fieldtrip/fieldtrip-20170514/fieldtrip-20170514/

%addpath misc/
addpath('C:\Chesney\EEG_data\From John');
addpath('thirdparty/fieldtrip/fieldtrip-20170514/fieldtrip-20170514/');

filelist = {'OA15','OA18','OA21','OA22','OA23','OA24'};

for fileno = (1:length(filelist))
    
    % Saved EEG structure after loading into EEGLAB and using the channel lookup function (gui)
    filestem = filelist{fileno};
    EEG = pop_loadset([filestem '_walk_long.set']);
    %EEG.chanlocs=pop_chanedit(EEG.chanlocs, 'load',{ '/Users/jsb/Desktop/EEG data/thirdparty/eeglab/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp', 'filetype', 'autodetect'});
    trigger = 3;
    EEG=pop_epoch(EEG,{trigger},[-1 1]);
    
    % Convert to Fieldtrip
    data = eeglab2fieldtrip( EEG, 'preprocessing', 'none' );

    % Sort out channel labels
    %data.elec.label{strcmp(data.label,'BIP3')} = 'acc1';
    data.label{strcmp(data.label,'BIP3') | strcmp(data.label,'ACC1')} = 'acc1';
    %data.elec.label{strcmp(data.label,'BIP4')} = 'acc2';
    data.label{strcmp(data.label,'BIP4') | strcmp(data.label,'ACC2')} = 'acc2';
    %data.elec.label{strcmp(data.label,'BIP5')} = 'acc3';
    data.label{strcmp(data.label,'BIP5') | strcmp(data.label,'ACC3')} = 'acc3';
    
    % Chuck initial portion of datasets
    if ( 0 )
        kk = round(chuck_initial*data.fsample):size(data.trial{1},2);
        data.time{1} = data.time{1}(kk);
        data.trial{1} = data.trial{1}(:,kk);
    end

    % Initial look at data
    if ( 0 )
        stackplotsc(data.trial{1}');
    end

    % Continuous data pre-processing
    cfg = [];
    cfg.channel = 'all';
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    data = ft_preprocessing(cfg, data);

    % Epoch data
    if ( 0 )
        cfg = [];
        triallen = 2.0;
        trialN = ceil(triallen*data.fsample);
        trialsamples = (1:trialN:length(data.time{1}))';
        cfg.trl = [ trialsamples(1:end-1) trialsamples(2:end)-1 zeros(length(trialsamples)-1,1) ];
        data = ft_redefinetrial(cfg, data);
    end

    % Visual artifact rejection (EEG channels only)
    cfg=[];     % Keep non-eeg channels
    cfg.channel = {'eog','acc1','acc2','acc3','eegref'};
    dataOther = ft_selectdata(cfg,data);
    cfg=[];
    cfg.channel = {'all','-eog','-acc1','-acc2','-acc3','-eegref'};
    dataEEG = ft_selectdata(cfg,data);
    cfg=[];
    cfg.keepchannel = 'repair';
    cfg.neighbours = '';
    dataEEG = ft_rejectvisual( cfg,  dataEEG );
    cfg=[];     % Merge datasets
    data = ft_appenddata(cfg,dataEEG,dataOther);
    
    %% Regress out Acc

    d = cat(3,data.trial{:});

    single_trial = true;
    dims = size(d);
    if ( single_trial )
        d = reshape( permute(d,[2 1 3]), [size(d,2) size(d,1)*size(d,3)] );
    end

    f = d;
    for seg = (1:size(d,3))
        X=[ ones(size(d,2),1) d((end-2):end,:,seg)' ];
        for ch = (1:(size(d,1)-3))
            b = regress(d(ch,:,seg)',X);
            f(ch,:,seg) = d(ch,:,seg) - (X*b)';
        end
    end

    if ( single_trial )
        f = permute( reshape( f, [dims(2) dims(1) dims(3)] ), [2 1 3] );
    end

    for seg = (1:size(f,3))
        data.trial{seg} = f(:,:,seg);
    end

    %% ICA decomposition

    if ( 0 )

        cfg            = [];
        cfg.method     = 'runica';
        comp           = ft_componentanalysis(cfg, data);

        % IC topoplot

        cfg           = [];
        cfg.component = [1:45];       % specify the component(s) that should be plotted
        cfg.layout    = 'easycapM25.lay'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)

        % ICA browser

        cfg          = [];
        cfg.channel  = (1:45); % components to be plotted
        cfg.viewmode = 'component';
        cfg.layout   = 'easycapM25.lay'; % specify the layout file that should be used for plotting
        cfg2 = ft_databrowser(cfg, comp)

        % Reject components

        cfg = [];
        cfg.component = [];%[ 22 29 33 38 40 42 44 ];
        data0=data;
        data = ft_rejectcomponent(cfg, comp, data);

    end

    %% Downsample

    cfg = [];
    cfg.resamplefs  = 128;
    cfg.detrend     = 'yes';
    data = ft_resampledata(cfg, data);

    save([filestem '_walk_long_epoch' num2str(trigger) '_ds'],'data');

    %% Rereference

    channel = {'all','-eog','-acc1','-acc2','-acc3','-eegref'};
    chanlist = ft_channelselection(channel,data);

    kklist = ismember(data.label,chanlist);
    for seg = (1:length(chanlist))
        data.trial{seg} = data.trial{seg} - mean(data.trial{seg}(kklist,:),1);
    end

    %% Look at time-frequency spectra

    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:64;                         % analysis 2 to 30 Hz in steps of 2 Hz 
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.2;   % length of time window = 0.2 sec
    cfg.toi          = -1.0:0.01:1.0;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
    TFRhann = ft_freqanalysis(cfg, data);

    if ( 0 )
        cfg = [];
        cfg.layout   = 'easycapM25.lay';
        cfg.ylim = [13 30];
        ft_topoplotTFR( cfg, TFRhann );
    end
    
    save([filestem '_walk_long_epoch' num2str(trigger) '_tfr'],'TFRhann');

end

return;

%%

cfg = [];
cfg.layout   = 'easycapM25.lay';
cfg.ylim = [13 30];

trigger=[];
trans = @(x) repmat( log(nanmax(x,[],3)-nanmin(x,[],3)), [1 1 size(x,3)] );

clear('TFR'); tfrmean = [];
for fileno = (1:length(filelist))
    
    % Saved EEG structure after loading into EEGLAB and using the channel lookup function (gui)
    filestem = filelist{fileno};
    TFR(fileno)=load([filestem '_walk_long_epoch' num2str(trigger) '_tfr'],'TFRhann');
    %TFR(fileno).TFRhann.powspctrm = trans(TFR(fileno).TFRhann.powspctrm);
    
    subplot(3,3,fileno);
    ft_topoplotTFR( cfg, TFR(fileno).TFRhann );
    
    tfrmean(:,:,:,fileno) = TFR(fileno).TFRhann.powspctrm;
    title(filelist(fileno));
end

subplot(3,3,8);
TFRhann = TFR(1).TFRhann;
TFRhann.powspctrm = nanmean(log10(tfrmean),4);
ft_topoplotTFR( cfg, TFRhann );
title('Average');

subplot(3,3,9);
TFRhann = TFR(1).TFRhann;
TFRhann.powspctrm = nanmean(tfrmean,4)./(nanstd(tfrmean,[],4)/sqrt(size(tfrmean,4)));
ft_topoplotTFR( cfg, TFRhann );
title('t-statistic');

%%

metric = TFRhann;
metric.powspctrm = repmat(log(max(metric.powspctrm,[],3)./min(metric.powspctrm,[],3)),[1 1 size(metric.powspctrm,3)]);

cfg = [];
cfg.layout   = 'easycapM25.lay';
cfg.ylim = [13 30];
ft_topoplotTFR( cfg, metric );

%% Plot power spectra

cfg = [];
cfg.baseline = [-.5 .5];
cfg.baselinetype = 'relative';
cfg.xlim = [-0.5 0.5];
cfg.zlim = [0.65 3];
cfg.layout   = 'easycapM25.lay';
cfg.showlabels = 'yes';

figure;
ft_multiplotTFR(cfg, TFRhann);

%% Plot power spectra

cfg = [];
cfg.baseline = [-1 1];
cfg.baselinetype = 'relative';
cfg.xlim = [-0.4 0.4];

cfg.channel = 'acc1';
cfg.showlabels = 'yes';
figure; ft_singleplotTFR(cfg, TFRhann);
cfg.channel = 'acc2';
figure; ft_singleplotTFR(cfg, TFRhann);
cfg.channel = 'acc3';
figure; ft_singleplotTFR(cfg, TFRhann);

cfg.channel = 'C3';
figure; ft_singleplotTFR(cfg, TFRhann);

%%

cfg = [];
cfg.baseline     = [];
cfg.baselinetype = 'absolute';  
cfg.maskstyle    = 'saturation';	
cfg.channel      = 'Cz';
figure 
ft_singleplotTFR(cfg, TFRhann);

cfg = [];
cfg.baseline     = [-0.5 -0.1];	
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.9 1.3];   
cfg.zlim         = [-1.5e-27 1.5e-27];
cfg.ylim         = [15 20];
cfg.marker       = 'on';
figure 
ft_topoplotTFR(cfg, TFRhann);

%% Freq_Resp_DAQ_TACS_Control
% V.3. - ERhodes 29/03/23
% This script will connect to the National Instruments USB 6216 DAQ, which
% in turn will be connected to the input of the Soterix Medical Stimulator we will be
% using to deliver retinal stimulation.
% The script will generate the waveform to be sent to the stimulator, the
% participant will then respond whether they saw any phosphenes or not. The
% threshold detection will then be aided using binary search.
%% Pre-amble
close all
clearvars -except s
clc

%% Directories
cF = pwd; % Working directory
% Functions
addpath(fullfile(cF,'functions'))
% Already created and populated with subject folders which will also include
% the condition order per subject:
dF = fullfile(cF, '0 - raw_data');
% Load Subject Order
answer = inputdlg({'Subject ID: ', 'Site No (1 or 2): ','Freq No: ' } ,'Subject to Test',[1,35],{'99','1','1'});

sbj_num = str2double(answer{1});
site_2_strt = str2double(answer{2});
freq_2_strt = str2double(answer{3});
sbj_id = num2str(sbj_num);
sdF = fullfile(dF,['s' sbj_id]);
load(fullfile(sdF,['s' sbj_id '_cond_order.mat']))

%% Settings
%%%%%% DAQ Settings %%%%%%%
Fs = 10e3; %10kHz

%%%% Stimulation Waveform %%%%
beat_F = 16; %16Hz is the optimal phosphene frequency
% 1s Per Amp Step, 30s total duration
amp_step = 0.25/30;
min_amps = [0,0.25,0.5,0.75];
max_amps = [0.25,0.5,0.75,1];
num_ramps = numel(min_amps);
for i = num_ramps:-1:1
    tmp_amp_ind = min_amps(i)+amp_step:amp_step:max_amps(i);
    tmp_amp_arr = ones(1,30*Fs);
    for i2 = 30:-1:1
        tmp_amp_arr(((i2-1)*Fs)+1:i2*Fs) = tmp_amp_ind(i2);
    end
    ramp_amps(i,:) = tmp_amp_arr;
    clear tmp_amp_ind tmp_amp_arr
end

%%% Binary Search Waveform Is Longer %%%
ramp_up_dur = 3; %s
stim_dur = 5; %s
ramp_down_dur = 2; %s
total_dur = ramp_up_dur+stim_dur+ramp_down_dur;
static_amp = [linspace(0,1,ramp_up_dur*Fs),ones(1,stim_dur*Fs),linspace(1,0,ramp_down_dur*Fs)];

%% Connect to DAQ
%New
if ~exist('s','VAR')
    s = daq('ni');
    out_ch = addoutput(s,'Dev1','ao0','Voltage');
    in_ch = addinput(s,'Dev1','ai7','Voltage');
    s.Rate = Fs;
end
flush(s);
%% Display Subject Info
clc
disp(['%%% Subject ID : ' num2str(sbj_id) ' %%%'])

%% Start Stimulation Experiment
pain_det =[];
for stim = site_2_strt:2
    if stim == 1
        disp(['First stimulation site will be the ' sbj_site_ord{1}])
    else
        disp(['Change stimulation site to the ' sbj_site_ord{2}])
    end
    pause(20)
    questdlg('EXPERIMENTER - Ready to Begin?',sbj_site_ord{stim},...
        'Yes','Yes');
    flush(s)
    if strcmp(sbj_site_ord{stim}, 'arm')
        carr_f_ord = sbj_peripheral_freq_ord;
        quest_prompt = 'Did you feel any sensation on your arm?';
    elseif strcmp(sbj_site_ord{stim}, 'retina')
        carr_f_ord = sbj_phosphene_freq_ord;
        quest_prompt = 'Do you see any flickering light?';
    end
    
    for freq = freq_2_strt:numel(carr_f_ord)
        disp(['%% Frequency No: ' num2str(freq) ' %%'])
        carr_F = carr_f_ord(freq);
        %%%% Start With Ramping Stim Amp %%%%
        if carr_F == 0
            out_wave = ISN_GenSine(beat_F,30,Fs);
        else
            out_wave = ISN_GenAM(beat_F,carr_F,30,Fs);
        end
        user_stop = 0;
        for i = 1:num_ramps
            disp(['% Ramping stim ' num2str(i) '/4 %'])
            stim_out = out_wave.*ramp_amps(i,:);
            stim_out = [zeros(1,5*Fs), stim_out];
            preload(s,stim_out');
            if i > 1
                f = msgbox('Pause for 10s');
                pause(8)
                delete(f)
            end
            str_qst = questdlg('When you are ready, press the ENTER key','Start Stimulation',...
                'Start','Cancel','Start');
            switch str_qst
                case 'Cancel'
                    error('User stopped experiment')
                    flush(s)
            end
            start(s)
            tic
            pause(2)
            button = questdlgtimeout(33, quest_prompt,'Stimulation Detection',...
                'Yes','Yes');
            but_stop = toc;
            stop(s)
            
            if but_stop < 35
                user_stop =1;
            end
            
            if user_stop == 1
                tmp_dat = read(s,'all');
                flush(s);
                stop_pnt = floor(size(tmp_dat,1)/Fs)-5;
                init_thresh = ramp_amps(i,(stop_pnt+.5)*Fs);
                if strcmp(sbj_site_ord{stim}, 'retina')
                    pain_det = [];
                    phos_qst = questdlg('Participant stopped due to miscomfort or seeing phosphene?','Phosphene Check',...
                        'Phosphene','Uncomfortable','Phosphene');
                    switch phos_qst
                        case 'Uncomfortable'
                            user_stop = 0;
                            pain_det =init_thresh;
                    end
                end
                
                break;
            else
                flush(s);
            end
            
            
            
            if i == num_ramps
                uns_qst = questdlg('Was participant unsure at the end of block?','Participant Check',...
                    'Yes','No','No');
                switch uns_qst
                    case 'Yes'
                        user_stop = 1;
                        init_thresh = 1;
                end
            end
        end
        
        
        %%%%% If a rough threshold is found, binary search %%%%%
        if user_stop == 1
            disp('%%% Starting Binary Search %%%%%')
            yes_count = 0;
            out_wave = out_wave(1:size(static_amp,2));
            % Binary Search
            if init_thresh == 1
                min_amp = 1-.125;
                max_amp = 1+.125;
            else
                min_amp = .5*init_thresh;
                max_amp = 1.5*init_thresh;
            end
            strt_amp = max_amp/2;
            cur_amp = strt_amp;
            lst_det = max_amp;
            low_thresh = min_amp;
            diff_amp = lst_det - low_thresh;
            i = 0;
            while diff_amp > .005
                i = i+1;
                stim_out = out_wave*cur_amp;
                stim_out = stim_out.*static_amp;
                flush(s)
                preload(s,stim_out')
                str_qst = questdlg('When you are ready, press the ENTER key','Start Stimulation',...
                    'Start','Cancel','Start');
                switch str_qst
                    case 'Cancel'
                        error('User stopped experiment')
                        flush(s)
                end
                start(s)
                pause(10.5)
                stop(s)
                stim_det = questdlg(quest_prompt,'Stim Detection',...
                    'No','Yes','No');
                switch stim_det
                    case 'No'
                        low_thresh = cur_amp;
                        cur_amp = .5 * (cur_amp+lst_det);
                    case 'Yes'
                        lst_det = cur_amp;
                        cur_amp = .5 * (cur_amp+low_thresh);
                        yes_count = yes_count + 1;
                end
                diff_amp = lst_det - low_thresh;
            end
            flush(s)
            det_thresh = lst_det;
            if yes_count > 0
                save(fullfile(sdF,[sbj_site_ord{stim} '_' num2str(carr_F) 'Hz.mat']),'det_thresh')
            else
                error(['Binary search failed: Rerun Site :' num2str(stim) ', Freq :' num2str(freq)])
            end
        else
            % No Detection
            det_thresh = NaN;
            if ~isempty(pain_det)
                save(fullfile(sdF,[sbj_site_ord{stim} '_' num2str(carr_F) 'Hz.mat']),'det_thresh','pain_det')
            else
                save(fullfile(sdF,[sbj_site_ord{stim} '_' num2str(carr_F) 'Hz.mat']),'det_thresh')
            end
        end
        
    end
end

woo = msgbox('Experiment complete! Thank You!');

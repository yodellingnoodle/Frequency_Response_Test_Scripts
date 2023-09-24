%% Preamble
close all
clear
clc

%% Directories
cF = pwd;
addpath(fullfile(cF,'functions'));
dF = fullfile(cF,'0 - raw_data');

%% Conditions
% We will do peripheral (arm) and transorbital (retinal) stimulation
stim_site ={'arm','retina'};
num_site = numel(stim_site);
% A range of carrier frequencies will be used, beat frequency will be 16Hz
car_frex =[0,50,75,100,250,500,750,1000,2500,5000];
num_car_f = numel(car_frex);

%% Randomise Carrier Frequencies
[~,rand_freq_ord1] = ISN_ConditionBalancing(num_car_f,0);
[~,rand_freq_ord2] = ISN_ConditionBalancing(num_car_f,0);
peripheral_freq_ord = [rand_freq_ord1;rand_freq_ord2];
peripheral_freq_ord = car_frex(peripheral_freq_ord);
phosphene_freq_ord = [rand_freq_ord2;rand_freq_ord1];
phosphene_freq_ord = car_frex(phosphene_freq_ord);
num_sbjcts = size(peripheral_freq_ord,1);
clear rand_freq_ord1 rand_freq_ord2

%% Randomise Stim Site Order
tmp = ones(num_sbjcts,2);
tmp(1:num_sbjcts/2,2) = tmp(1:num_sbjcts/2,2) *2;
tmp((num_sbjcts/2)+1:end,1) = tmp((num_sbjcts/2)+1:end,1) *2;
site_ord = tmp(randperm(num_sbjcts),:);
site_ord = stim_site(site_ord);
clear tmp

%% Save All Subject Order
save(fullfile(dF,'all_sbj_condition_order.mat'),'site_ord',...
    'peripheral_freq_ord','phosphene_freq_ord')

%% Create Subject Folders
% Define SXX folder
sxF = fullfile(dF,'sxx','docs');
for n = num_sbjcts:-1:1
    % Make Main Subject Folder
    sF = fullfile(dF,['s' num2str(n)]);
    if ~exist(sF,'dir')
        mkdir(sF)
    end
    % Save Subject Condition Order
    sbj_peripheral_freq_ord = peripheral_freq_ord(n,:);
    sbj_phosphene_freq_ord = phosphene_freq_ord(n,:);
    sbj_site_ord = site_ord(n,:);
    save(fullfile(sF,['s' num2str(n) '_cond_order.mat']),'sbj_site_ord',...
        'sbj_peripheral_freq_ord','sbj_phosphene_freq_ord')
    clear sbj_peripheral_freq_ord sbj_phosphene_freq_ord sbj_site_ord
    
    % Copy Documents from SXX to Current Subject
    doc_fld = fullfile(sF,'docs');
    mkdir(doc_fld)
    copyfile(fullfile(sxF,'s_screening.docx'),fullfile(doc_fld,['s' num2str(n) '_screening.docx']))
    copyfile(fullfile(sxF,'Participant_Information_Sheet_v1.4_18_11_2020.docx'),...
        fullfile(doc_fld,'Participant_Information_Sheet_v1.4_18_11_2020.docx'))
    copyfile(fullfile(sxF,'Consent_Form.pdf'),fullfile(doc_fld,'Consent_Form.pdf'))
end

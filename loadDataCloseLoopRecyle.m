%%% Used to analyze Place cell on LS data
%%% This pipeline will define trials
%%% 202210_HJY
%%% data location
%%% cell_resp: /10.10.49.10/hjy/Light Sheet Data/20221018_1_2_8f_5dpf/20221018_1_2_8f_5dpf_20221018_205617/registered/cell_resp.stackf
%%% ephys: /10.10.49.10/hjy/Light Sheet Data/20221018_1_2_8f_5dpf/20221018_1_2_8f_5dpf_20221018_205617/ephys/20221018_1_2_8f_5dpf.10chFlt
%%%  Used for VR version
clear all;close all;

%% Ephy data loading
% ephy_dir = 'Z:\hjy\Light Sheet Data\20221029_3_1_8f_6dpf\20221029_3_1_8f_6dpf_20221030_001013\ephys';

ephy_dir = 'Z:\hjy\Light Sheet Data\20221029_2_1_8f_6dpf\20221029_2_1_8f_6dpf_20221029_204707\ephys';


%ephy_data = ephy_read_data_13(fullfile(ephy_dir,'20221018_1_2_8f_5dpf.10chFlt'));
ephy_data = read_data(fullfile(ephy_dir,'20221029_2_1_8f_6dpf.10chFlt')); %20221018_1_2_8f_5dpf

% useless = {'ch4','fishx','stimpara','fishy'}; 
useless = {'ch4','stimpara'}; 
ephy_data = rmfield(ephy_data,useless);

%posY-0.201/figure -0.201
% 0.4078 fish head
% 
%% Light Sheet data loading
disp('Loading LS data......');
LSinput_dir = '/Volumes/10.10.49.10/hjy/Light Sheet Data/20221029_3_1_8f_6dpf/20221029_3_1_8f_6dpf_20221030_001013\registered'
load(fullfile(LSinput_dir,'cell_resp_dim_processed.mat'));   % cell_resp_dim
cell_resp = read_LSstack_fast_float(fullfile(LSinput_dir,'cell_resp_processed.stackf'),cell_resp_dim);
load(fullfile(LSinput_dir,'cell_info_processed.mat'));
%cell_resp_raw = read_LSstack_fast_float(fullfile(LSinput_dir,'cell_resp.stackf'),cell_resp_dim); % cell_resp_raw

% %% data preprocessing
% filtdata_l = filter_data(ephy_data.ch2); % channel 2 is left ephy
% bouts = extract_bouts(filtdata_l);
%% align ephy data with light-sheet data
% camera_trigger_t = find(ephy_data.ch3(1:end-1)<3.6 & ephy_data.ch3(2:end)>3.6)+1; % camera trigger in ephys
% frame_info = load(fullfile(LSinput_dir,'frame_info.txt'));  % processed frame
% camera_trigger_t = camera_trigger_t(frame_info(1):frame_info(2));
% 
% ephy_data.trigger = zeros(length(ephy_data.ch2),1); % ephys t - camera t Align with the camera trigger behind
% ephy_data.trigger(1:camera_trigger_t(1)-min(diff(camera_trigger_t))-1) = 0;
% ephy_data.trigger(camera_trigger_t(1)-min(diff(camera_trigger_t)):camera_trigger_t(1)) = 1;
% for i = 1:length(camera_trigger_t)-1
%     ephy_data.trigger(camera_trigger_t(i)+1:camera_trigger_t(i+1)) = i+1;
% end

% plane_t = find(diff(ephy_data.ch3)>2)+1;
% xx = diff(plane_t);
% delt = find(xx>300)+1;
% camera_out = plane_t(delt(2:end));
% 
% ephy_data.trigger = zeros(length(ephy_data.ch2),1); % ephys t - camera t Align with the camera trigger behind
% ephy_data.trigger(1:camera_out(1)-min(diff(camera_out))-1) = 0;
% ephy_data.trigger(camera_out(1)-min(diff(camera_out)):camera_out(1)) = 1;
% for i = 1:length(camera_out)-1
%     ephy_data.trigger(camera_out(i)+1:camera_out(i+1)) = i+1;
% end

%% align ephy data with light-sheet data
plane_t = find(diff(ephy_data.ch3)>2)+1;
xx = diff(plane_t);
delt = find(xx>300)+1;
camera_trigger_t = plane_t(delt(2:end));

ephy_data.trigger = zeros(length(ephy_data.ch2),1); % ephys t - camera t Align with the camera trigger behind
ephy_data.trigger(1:camera_trigger_t(1)-min(diff(camera_trigger_t))-1) = 0;
ephy_data.trigger(camera_trigger_t(1)-min(diff(camera_trigger_t)):camera_trigger_t(1)) = 1;
for i = 1:length(camera_trigger_t)-1
    ephy_data.trigger(camera_trigger_t(i)+1:camera_trigger_t(i+1)) = i+1;
end
%% Trial defing
trial_length = 1;
ephy_data.trialY = ephy_data.posY;% +0.5+0.072/2;%0.201-0.072; %% calibrate the screen with the FOV of the chamber
ephy_data.trialY = mod(ephy_data.trialY,trial_length);
%% align position/context
posi = [];
context = [];
%[ ) forward mean 
% ephy_data.trialY
for i =1:size(cell_resp,2)
    indePosi0 = camera_trigger_t(i);
%     indePosi1 = camera_trigger_t(i+1);
    posi = [posi nanmean(ephy_data.trialY(indePosi0))];
    context = [context ephy_data.landm(indePosi0)]; 
end
%%
%correct landmark
contCorr = context;
for i =2:21430
    if posi(i)>=posi(i-1)
        contCorr(i) = contCorr(i-1);        
    end
end
% ylim([-0.1,1.1])
%% combine cell response with trial information
cellRespLand = [cell_resp;posi;contCorr];
%% save data
save('cellRespLand20221029_3_1.mat',...
    'cell_resp',...
    'context',...
    'posi',... 
    '-v7.3')














%% bin data by the fish's location
% By Guangyao Qi 2022-11-02

%% load data
clc
clear
load('cellRespLand20221029_3_1.mat')
%% test temporal feature
trialInde = [1];
t = 1;
for n =2:length(context)
    if context(n) >= context(n-1)
        trialInde = [trialInde,t];
    else
        t = t + 1;
        trialInde = [trialInde,t];
    end
end
%% 
delaPeri = [];
cond0Peri = [];
cond1Peri = [];
cond2Peri = [];
for i =1:length(unique(trialInde))
    delaPeri = [delaPeri, sum(trialInde==i & context == 3)];
    cond0Peri = [cond0Peri, sum(trialInde==i & context == 0)];
    cond1Peri = [cond1Peri, sum(trialInde==i & context == 1)];
    cond2Peri = [cond2Peri, sum(trialInde==i & context == 2)];
end
% hist()
badTria = [1,20,55];%3_1
% badTria = [1,18,19,26,54];%2_1
cond0Peri(badTria)= [];
cond1Peri(badTria) = [];
cond2Peri(badTria) = [];
%% bin by 5/10/15
% parpool(8)
tic
binsNumb = 5;
cell_resp_bins = nan(size(cell_resp,1),binsNumb*length(unique(trialInde))*3);
for n =1:size(cell_resp,1)
    n
    for t =1:length(unique(trialInde))
        for i =1:binsNumb
            for j =0:2
                cell_resp_bins(n,i+j*binsNumb+(t-1)*binsNumb*3) = nanmean(cell_resp(n,...
                    trialInde==t &...
                    context == j & ...
                    posi>=(i-1)*1/binsNumb &...
                    posi<i*1/binsNumb));
            end
        end
    end
end
toc
%% save bin data
save('cellRespLand20221029_3_1_5Bins.mat',...
    'cell_resp_bins')







































function pwr_95CI = calcRWABootstrap(pwr,trialtype)
%pwr is time x freq x trial
%trialtype is 1 x ntrials or scalar

ntime = size(pwr,1);
nfreq = size(pwr,2);
ntrials = size(pwr,3);
trials1 = find(trialtype==-1);
ntrials1 = length(trials1);
trials2 = find(trialtype==1);
ntrials2 = length(trials2);
B = 1000; %number of bootstrap iterations

if isscalar(trialtype)
    boottype = 1; %pwr
else
    boottype = 0; %ttest
end

rng('shuffle'); %seed the random stream with clock time
pwr_boot = nan(ntime,nfreq,B);
parfor k=1:B
    if boottype
        idx = randi(ntrials, [1 ntrials]); %randomly sample with replacement
        pwr_boot(:,:,k) = mean(pwr(:,:,idx),3,'omitnan'); %mean across trials
    else
        idx1 = trials1(randi(ntrials1, [1 ntrials1])); %randomly sample with replacement
        idx2 = trials2(randi(ntrials2, [1 ntrials2])); %randomly sample with replacement
        tnum = squeeze(mean(pwr(:,:,idx1),3,"omitnan")-mean(pwr(:,:,idx2),3,"omitnan")); %trial avg diff (time x freq)
        tdenom = sqrt(std(pwr(:,:,idx1),0,3,"omitnan").^2./ntrials1+std(pwr(:,:,idx2),0,3,"omitnan").^2./ntrials2);
        pwr_boot(:,:,k) = tnum./tdenom; %time x freq;
    end
end

pwr_95CI = prctile(pwr_boot,[2.5,97.5],3);
function PM = calcRWAITPCDiffPerm(cfs,pwr,trialtype,nperm)

if isempty(cfs)
    ntime = size(pwr,1);
    nfreq = size(pwr,2);
else
    ntime = size(cfs,1);
    nfreq = size(cfs,2);
end
ntrials1 = sum(trialtype==-1);
ntrials2 = sum(trialtype==1);

% rng('shuffle'); %seed the random stream with clock time
rng(1, 'twister'); %should use this to get reproducible distributions 
PM = nan(ntime,nfreq,nperm);
% wait_msg = parfor_wait(nperm);
parfor k=1:nperm
    % wait_msg.Send;
    ttshuf = trialtype(randperm(ntrials1+ntrials2));
    if isempty(cfs)
        tnum = squeeze(mean(pwr(:,:,ttshuf==-1),3)-mean(pwr(:,:,ttshuf==1),3)); %trial avg diff (time x freq)
        tdenom = sqrt(std(pwr(:,:,ttshuf==-1),0,3).^2./ntrials1+std(pwr(:,:,ttshuf==1),0,3).^2./ntrials2);
        PM(:,:,k) = tnum./tdenom;
    else
        ppc1 = (abs(sum(cfs(:,:,ttshuf==-1),3)).^2-ntrials1)./(ntrials1.*(ntrials1-1));
        ppc2 = (abs(sum(cfs(:,:,ttshuf==1),3)).^2-ntrials2)./(ntrials2.*(ntrials2-1));
        PM(:,:,k) = ppc1-ppc2;
    end
end
% wait_msg.Destroy;
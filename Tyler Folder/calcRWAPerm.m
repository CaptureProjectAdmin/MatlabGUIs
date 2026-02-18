function PM = calcRWAPerm(pwr,bpwr,trialtype,pwr_z)
%pwr (time x freq x trial), bpwr (1 x freq) -> baseline power, trialtype is
%-1 or 1 for calculating difference between specgrams or 0 for a single
%specgram. pwr_z is mean/std for applying zscore (2 x freq x trial).

ntime = size(pwr,1);
nfreq = size(pwr,2);
ntrials = size(pwr,3);
ntrials1 = sum(trialtype==-1);
ntrials2 = sum(trialtype==1);
nperm = 1000;

if isscalar(trialtype)
    permtype = 1; %shift
else
    permtype = 0; %ttest
end

%This should speed up the permutations by skipping the zscore calculation
%when the values are 0/1.
if all(pwr_z(1,:)==0) && all(pwr_z(2,:)==1)
    zflag = false;
else
    zflag = true;
end

rng('shuffle'); %seed the random stream with clock time
PM = nan(ntime,nfreq,nperm);
% wait_msg = parfor_wait(nperm);
parfor m=1:nperm %permutations
%     wait_msg.Send;
    if permtype %shift
        cutpoint = randi(ntime,[1,ntrials]);
        pwr_shift = zeros(size(pwr));
        for k=1:ntrials
            pwr_shift(:,:,k) = circshift(pwr(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
        end
        npwr = 10*log10(pwr_shift./bpwr);
        if zflag
            npwr = (npwr-pwr_z(1,:,:))./pwr_z(2,:,:);
        end
        PM(:,:,m) = squeeze(mean(npwr,3,"omitnan")); %normalized trial avg power in dB (take trial avg in log space)
    else %ttest
        ttshuf = trialtype(randperm(ntrials));
        tnum = squeeze(mean(pwr(:,:,ttshuf==-1),3,"omitnan")-mean(pwr(:,:,ttshuf==1),3,"omitnan")); %trial avg diff (time x freq)
        tdenom = sqrt(std(pwr(:,:,ttshuf==-1),0,3,"omitnan").^2./ntrials1+std(pwr(:,:,ttshuf==1),0,3,"omitnan").^2./ntrials2);
        PM(:,:,m) = tnum./tdenom; %time x freq;
    end
end
% wait_msg.Destroy;
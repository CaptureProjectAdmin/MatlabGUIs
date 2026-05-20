function PM = calcRWAITPCPermPwr(cfs)

ntime = size(cfs,1);
nfreq = size(cfs,2);
ntrials = size(cfs,3);
nperm = 1000;
rng('shuffle'); %seed the random stream with clock time
PM = nan(ntime,nfreq,nperm);
% wait_msg = parfor_wait(nperm);
parfor m=1:nperm %permutations
    % wait_msg.Send;
    cutpoint = randi(ntime,[1,ntrials]);
    cfs_shift = zeros(size(cfs));
    for k=1:ntrials
        cfs_shift(:,:,k) = circshift(cfs(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
    end
    PM(:,:,m) = mean(cfs_shift,3,'omitnan'); %mean across trials
end
% wait_msg.Destroy;
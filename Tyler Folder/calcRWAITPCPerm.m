function PM = calcRWAITPCPerm(cfs)

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
    cfs_shift = complex(zeros(size(cfs)));
    for k=1:ntrials
        cfs_shift(:,:,k) = circshift(cfs(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
    end
    % itpc = exp(1i*(angle(cfs_shift)));
    itpc = cfs_shift./abs(cfs_shift); %same as above but faster
    itpc = mean(itpc,3,'omitnan'); %mean across trials in complex
    % itpc = smoothdata(itpc,1,'movmean',3); %smoothing across time in complex
    PM(:,:,m) = abs(itpc);
end
% wait_msg.Destroy;
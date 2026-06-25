function PM = calcRWAITPCPerm(cfs,pwr,nperm)

if isempty(cfs)
    ntime = size(pwr,1);
    nfreq = size(pwr,2);
    ntrials = size(pwr,3);
else
    ntime = size(cfs,1);
    nfreq = size(cfs,2);
    ntrials = size(cfs,3);
end

% rng('shuffle'); %seed the random stream with clock time
rng(1, 'twister');
PM = nan(ntime,nfreq,nperm);
% wait_msg = parfor_wait(nperm);
parfor m=1:nperm %permutations
    % wait_msg.Send;
    cutpoint = randi(ntime,[1,ntrials]);
    if isempty(cfs)
        pwr_shift = zeros(size(pwr));
        for k=1:ntrials
            pwr_shift(:,:,k) = circshift(pwr(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
        end
        PM(:,:,m) = mean(pwr_shift,3); %mean across trials (10*log10 already applied)
    else
        cfs_shift = complex(zeros(size(cfs)));
        for k=1:ntrials
            cfs_shift(:,:,k) = circshift(cfs(:,:,k),cutpoint(k),1); %time x freq x trial (shift time)
        end
        PM(:,:,m) = abs(mean(cfs_shift,3)); %mean across trials in complex
    end
end
% wait_msg.Destroy;
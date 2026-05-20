function itpc_band_boot = calcRWAITPCBootstrap(itpc_cfs)

rng('shuffle'); %seed the random stream with clock time
N = size(itpc_cfs,3); %number of trials
B = 1000; %number of bootstrap iterations
itpc_band_boot = nan(size(itpc_cfs,1), B);
parfor k=1:B
    idx = randi(N, [1 N]); %randomly sample with replacement
    itpc_tf_b = abs(mean(itpc_cfs(:,:,idx), 3, 'omitnan')); %mean across trials in complex
    itpc_band_boot(:, k) = mean(itpc_tf_b, 2, 'omitnan'); %mean across freq in real
end
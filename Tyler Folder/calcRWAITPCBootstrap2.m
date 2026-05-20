function [itpc_band,itpc_band_boot] = calcRWAITPCBootstrap2(itpc_cfs)

% ---- Inputs ----
% itpc_cfs: time x freq x trial (unit phasors)
% B:       number of bootstrap iterations

B = 1000;

[T, Fband, N] = size(itpc_cfs); % time x Fband x trial

% Reshape so trials are rows (observations)
X = permute(itpc_cfs, [3 1 2]);                % trial x time x Fband
X = reshape(X, N, T*Fband);                  % trial x (time*Fband)

% Statistic function: returns itpc_band (time x 1) for a resampled set of trials
statfun = @(Xin) local_itpc_band_stat(Xin, T, Fband);

% Use parallel workers if you like
opts = statset('UseParallel', true);

% BCa 95% CI
itpc_band_boot = bootci(B, {statfun, X}, 'type','bca', 'alpha',0.05, 'Options', opts);
% itpc_band_boot is 2 x T : row 1 = lower, row 2 = upper

% Point estimate (computed the same way, using all trials)
itpc_band = statfun(X);

    % ---- helper ----
    function itpc_band = local_itpc_band_stat(Xin, T, Fband)
        % Xin: nTrialsSampled x (T*Fband)
        n = size(Xin,1);
        A = reshape(Xin.', T, Fband, n);         % time x Fband x n
        itpc_tf = abs(mean(A, 3, 'omitnan'));    % time x Fband
        itpc_band = mean(itpc_tf, 2, 'omitnan'); % time x 1
    end

end
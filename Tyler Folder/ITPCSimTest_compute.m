function S = ITPCSimTest_compute(sig, Fs, varargin)
%ITPCSimTest_compute  ITPC + power using the same pipeline as RWAnalysis2.
%
% Matches:
%   RWAnalysis2.updateMultTableEye  (Morse CWT, unit vectors, smooth, DS)
%   RWAnalysis2.plotMultTransITPCEye (ITPC = abs(mean(unit,3)), z-score perms)
%
% Uses existing helpers in this folder without modifying them:
%   morseSpecGram, calcRWAITPCPerm (or calcRWAITPCPerm_mex if present)

p = inputParser;
addParameter(p,'FPass',[2,120]);
addParameter(p,'SmoothSec',0.03);         % same as updateMultTableEye
addParameter(p,'DownsampleFactor',10);    % 250 -> 25 Hz
addParameter(p,'MaxFreqHz',32.5);         % same cutoff as plotMultTransITPCEye
addParameter(p,'NPerm',200);
addParameter(p,'PVal',0.05);
addParameter(p,'CorrectionType','pixel'); % 'pixel' or 'fdr'
addParameter(p,'FullEpochNorm',false);    % if true, 2nd median-norm (when ~fullwalknorm in GUI)
addParameter(p,'t',[]);                   % optional time vector (samples x 1)
parse(p,varargin{:});
opt = p.Results;

if ~isa(sig,'double'); sig = double(sig); end
[nsamp,ntrials] = size(sig);
if isempty(opt.t)
    t = ((0:nsamp-1)' ./ Fs);
    t = t - mean([t(1),t(end)]);          % center if caller did not pass t
else
    t = opt.t(:);
end

% ---- Morse spectrogram (same as updateMultTableEye) ----
% morseSpecGram returns time x freq x trial
[f,~,cfs] = morseSpecGram(sig,Fs,opt.FPass);
if ismatrix(cfs)
    cfs = reshape(cfs,size(cfs,1),size(cfs,2),1);
end

% Power path
pwr = abs(cfs).^2;
pwr = smoothdata(pwr,1,'movmean',Fs*opt.SmoothSec);
pwr = pwr ./ median(pwr,1,'omitnan');           % "full walk" / full-epoch median norm
pwr = pwr(1:opt.DownsampleFactor:end,:,:);

% ITPC unit-vector path
itpc_uv = cfs ./ abs(cfs);
itpc_uv = smoothdata(itpc_uv,1,'movmean',Fs*opt.SmoothSec);
itpc_uv = itpc_uv(1:opt.DownsampleFactor:end,:,:);

t_ds = t(1:opt.DownsampleFactor:end);

% Restrict to <= MaxFreqHz (same as plotMultTransITPCEye)
fkeep = f <= opt.MaxFreqHz;
f = f(fkeep);
itpc_cfs = itpc_uv(:,fkeep,:);
pwr_cfs = pwr(:,fkeep,:);

% ITPC
itpc = abs(mean(itpc_cfs,3));                   % abs of complex trial mean

% Power (dB)
if opt.FullEpochNorm
    pwr_cfs = pwr_cfs ./ median(pwr_cfs,1);     % epoch median norm (when ~fullwalknorm)
end
pwr_cfs = 10*log10(pwr_cfs + eps);
pwr_mean = mean(pwr_cfs,3);

% ---- Permutations (same as plotMultTransITPCEye) ----
permFun = localPermFun();
fprintf('ITPCSimTest: ITPC permutations (%d trials, %d perms)...\n',ntrials,opt.NPerm);
PM_itpc = permFun(itpc_cfs,[],opt.NPerm);
mPM = mean(PM_itpc,3);
sPM = std(PM_itpc,0,3);
zITPC = (itpc - mPM) ./ sPM;
zPM = (PM_itpc - mPM) ./ sPM;
zITPC_thresh = localThresh(zITPC,zPM,PM_itpc,itpc,opt);

fprintf('ITPCSimTest: power permutations (%d trials, %d perms)...\n',ntrials,opt.NPerm);
PM_pwr = permFun([],pwr_cfs,opt.NPerm);
mPMp = mean(PM_pwr,3);
sPMp = std(PM_pwr,0,3);
zPwr = (pwr_mean - mPMp) ./ sPMp;
zPMp = (PM_pwr - mPMp) ./ sPMp;
zPwr_thresh = localThresh(zPwr,zPMp,PM_pwr,pwr_mean,opt);

S = struct();
S.t = t;
S.t_ds = t_ds;
S.f = f;
S.sig = sig;
S.itpc_cfs = itpc_cfs;
S.itpc = itpc;
S.zITPC = zITPC;
S.zITPC_thresh = zITPC_thresh;
S.pwr_cfs = pwr_cfs;
S.pwr = pwr_mean;
S.zPwr = zPwr;
S.zPwr_thresh = zPwr_thresh;
S.ntrials = ntrials;
S.Fs = Fs;
S.nperm = opt.NPerm;
S.pval = opt.PVal;
S.correctiontype = opt.CorrectionType;
S.permtype = 'zscore';
end

% -------------------------------------------------------------------------
function fun = localPermFun()
if exist('calcRWAITPCPerm_mex','file') == 3
    fun = @calcRWAITPCPerm_mex;
else
    fun = @calcRWAITPCPerm;
end
end

function z_thresh = localThresh(zObs,zNull,PM,obs,opt)
z_thresh = zObs;
switch lower(opt.CorrectionType)
    case 'pixel'
        maxAbs = squeeze(max(max(abs(zNull),[],1),[],2));
        zcrit = prctile(maxAbs,100*(1-opt.PVal));
        z_thresh(abs(zObs) < zcrit) = 0;
    case 'fdr'
        pUpper = (1 + sum(PM >= obs,3)) ./ (opt.NPerm + 1);
        pLower = (1 + sum(PM <= obs,3)) ./ (opt.NPerm + 1);
        PV = min(1, 2 .* min(pUpper,pLower));
        QV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
        z_thresh(QV >= opt.PVal) = 0;
    otherwise
        error('Unknown CorrectionType: %s',opt.CorrectionType);
end
end

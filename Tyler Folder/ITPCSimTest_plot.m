function fH = ITPCSimTest_plot(S, meta, varargin)
%ITPCSimTest_plot  Match plotMultTransITPCEye layout: ITPC, power, ERP.
%
% Tile 1: z-scored ITPC contourf + white significance contour
% Tile 2: z-scored power contourf + white significance contour
% Tile 3: trial-average ERP +/- 95% SEM

p = inputParser;
addParameter(p,'CLim',[-10,10]);
addParameter(p,'YLimERP',[]);
parse(p,varargin{:});
opt = p.Results;

fH = figure('Position',[50,50,800,800],'Visible','on','Color','w');
tl = tiledlayout(3,1,'Parent',fH,'TileSpacing','compact','Padding','compact');

% ---- ITPC ----
aH = nexttile(tl,1);
colormap(aH,'jet');
hold(aH,'on');
contourf(S.t_ds,S.f,S.zITPC',100,'linecolor','none','parent',aH);
contour(S.t_ds,S.f,(S.zITPC_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
set(aH,'yscale','log','YTick',2.^(1:5),'yticklabel',2.^(1:5),'yminortick','off');
plot(aH,[0,0],[2,32],'--k','LineWidth',2);
axis(aH,[S.t_ds(1),S.t_ds(end),2,32]);
clim(aH,opt.CLim);
xlabel(aH,'sec');
ylabel(aH,'Hz');
ttl = sprintf(['ITPCSimTest theta (4-8 Hz), phase jitter \\pm%.0f deg, ' ...
    'f=%.1f\\pm%.1f Hz\\nn=%d, p<%.2f, %s, %s, onset@t=0 (abs %.0fs / %.0fs)'], ...
    rad2deg(meta.PhaseJitterRad), mean(meta.FreqHz), std(meta.FreqHz), ...
    S.ntrials, S.pval, S.permtype, S.correctiontype, meta.OnsetSec, meta.DurationSec);
title(aH,ttl);
cb = colorbar(aH);
cblims = [cb.Limits(1),0,cb.Limits(2)];
cb.Ticks = cblims;
cb.TickLabels = cblims;
ylabel(cb,'Zscore');

% ---- Power ----
aH = nexttile(tl,2);
colormap(aH,'jet');
hold(aH,'on');
contourf(S.t_ds,S.f,S.zPwr',100,'linecolor','none','parent',aH);
contour(S.t_ds,S.f,(S.zPwr_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
set(aH,'yscale','log','YTick',2.^(1:5),'yticklabel',2.^(1:5),'yminortick','off');
plot(aH,[0,0],[2,32],'--k','LineWidth',2);
axis(aH,[S.t_ds(1),S.t_ds(end),2,32]);
clim(aH,opt.CLim);
xlabel(aH,'sec');
ylabel(aH,'Hz');
cb = colorbar(aH);
cblims = [round(cb.Limits(1),1),0,round(cb.Limits(2),1)];
cb.Ticks = cblims;
cb.TickLabels = cblims;
ylabel(cb,'Zscore');
title(aH,'Power (z-scored vs circular-shift null)');

% ---- ERP ----
aH = nexttile(tl,3);
d = S.sig - mean(S.sig,1);
md = mean(d,2);
sd = std(d,0,2);
sd = sd ./ sqrt(size(S.sig,2)) * 1.96;
plot(aH,S.t,md);
hold(aH,'on');
plot(aH,S.t,md+sd,':k');
plot(aH,S.t,md-sd,':k');
plot(aH,[0,0],ylim(aH),'--k','LineWidth',1.5);
if ~isempty(opt.YLimERP)
    ylim(aH,opt.YLimERP);
end
ylabel(aH,'uV');
xlabel(aH,'sec');
title(aH,'Trial-average waveform (mean +/- 95% SEM)');
grid(aH,'on');

fH.UserData.S = S;
fH.UserData.meta = meta;
end

function dd = alignRWAITPCData(tsamp_np,ntp_np,d_np,ntp_gz)

dd = nan(length(tsamp_np),length(ntp_gz));
parfor m=1:length(ntp_gz)
    [~,idx] = min(abs(ntp_np-ntp_gz(m)));
    idx = idx+tsamp_np;
    if idx(1)>=1 && idx(end)<=length(d_np)
        dd(:,m) = d_np(idx);
    end
end
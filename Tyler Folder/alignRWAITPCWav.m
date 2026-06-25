function [dd,dd_it] = alignRWAITPCWav(tsamp_wv,ntp_wv,d_wv,ntp_gz,varargin)

d_it = complex([]);
if nargin>4
    d_it = varargin{1}; %itpc (complex, same size as d_wv)
end

dd_it = complex(nan(length(tsamp_wv),size(d_wv,2),length(ntp_gz)));
dd = nan(length(tsamp_wv),size(d_wv,2),length(ntp_gz));
parfor m=1:length(ntp_gz)
    [~,idx] = min(abs(ntp_wv-ntp_gz(m)));
    idx = idx+tsamp_wv;
    if idx(1)>=1 && idx(end)<=size(d_wv,1)
        dd(:,:,m) = d_wv(idx,:);
        if ~isempty(d_it)
            dd_it(:,:,m) = d_it(idx,:);
        end
    end
end
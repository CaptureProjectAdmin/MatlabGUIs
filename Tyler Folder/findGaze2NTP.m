function [ntp,fixD,sacD,Fs] = findGaze2NTP(Files,varargin)
%Open saccade/fixation files by Alireza and interpolates time for even
%sampling. This also produces separate data for saccades and fixations, which
%is different than the original implementation.

warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%get timestamps from csv file (need to adjust for drift)
drift_tbl = readtable(Files.drift_csv_file); %in units milliseconds
ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

gaze_tbl = readtable(Files.pupil_gaze_file,'Delimiter',',');

load(Files.gaze2_file,"sfdat");
sfdat.ntp_gaze.TimeZone = 'America/Denver';
ntp_alireza = datetime(sfdat.ntp_gaze,'timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
pt = str2double(regexp(Files.pt_dir,'\d$','match','once'));
wk = str2double(regexp(Files.walk_dir,'\d$','match','once'));
if any(pt==[3,5]) || ((pt==2) && (wk==5))
    ntp_alireza = ntp_alireza - hours(6);
end
d_fix_binary = sfdat.d_fix_binary;
d_sac_binary = sfdat.d_sac_binary;

ntp_orig = datetime(gaze_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');

if milliseconds(ntp_alireza(1)-ntp_orig(1))>1
    error('Original and Alireza timestamps do not match!');
end

ntp_orig = datenum(ntp_alireza)*60*60*24; %convert from days to seconds
ntp_orig = ntp_orig + ntp_offset;

Fs = 200; %there is some spiking in the sampling rate but mean is close to 200 (resampling to make consistent)
N = round((ntp_orig(end)-ntp_orig(1))*Fs); %number of samples based on ntp times and Fs
ntp = ((0:N-1)/Fs + ntp_orig(1))'; %time vector in sec that is evenly sampled at Fs and adjusted to match start of ntp

fixD = interp1(ntp_orig,double(d_fix_binary),ntp,'nearest','extrap'); %these are fixation counts so nearest is needed
sacD = interp1(ntp_orig,double(d_sac_binary),ntp,'nearest','extrap'); %these are fixation counts so nearest is needed








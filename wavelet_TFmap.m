function [Mag,taxis,faxis]=wavelet_TFmap(signal,num_of_columns,samplerate,f_low,f_high,f_step)
%%
% samplerate: in Hz
% fstep: frequency step for wavelet
% Tzu-Hao Harry Chao 2020/11/3
%%

data_range=1:length(signal); % data points of each time course
for i=1:num_of_columns
timecourse_data(:,i)=signal(data_range,i); % data matrix time courses in columns
forder=6;     % filter order
N=length(data_range);
spec = tfa_morlet(timecourse_data(:,i)', samplerate, f_low, f_high, f_step);
taxis=[1:N]/samplerate; % the x axis for plotting the time-freq figure
faxis=[f_low:f_step:f_high]; % the y axis for plotting the time-freq figure
Mag(:,:,i)=abs(spec);
end
% % %
% In the following steps, if the results are transformed by 10*log10,
% then the unit is dB/Hz, the original unit is "%/Hz" if the raw data
% had been transformed into percentage changes.
% % %
% plot the mean time-freq result (dB/Hz)
pcolor(taxis,faxis,(mean(Mag,3)))
shading flat;

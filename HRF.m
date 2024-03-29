function [HRF]=HRF(neuronal_activity,hemodynamics,sampling_rate,HRF_length)
%%
% calculate HRF via simultaneously recorded neuronal_activity and
% hemodynamics time-courses.
%
% neuronal_activity: neuronal activity time-course in a single column.
% hemodynamics: hemodynamics changes time-course in a single column.
% HRF_length: presumed HRF time length in seconds. Recommend starting with 25s.
%
% Example: HRF(GCaMP_timecourse,Rhodamine_timecourse,10,25,0)
% Tzu-Hao Harry Chao 2021/10/22
%%

hemodynamics=hemodynamics-min(hemodynamics);

HRF_length=HRF_length*sampling_rate+1; % sec to data points

clear X
raw_data_length=length(neuronal_activity);
X = zeros(raw_data_length,HRF_length);
temp = neuronal_activity(1:raw_data_length);
for i=1:HRF_length
X(:,i) = temp;
temp = [0;temp(1:end-1)];
end
X(:,HRF_length+1)=ones(raw_data_length,1)';
X(:,HRF_length+2)=linspace(0,1,raw_data_length)';

HRF=pinv(X(1:raw_data_length,:))*hemodynamics(1:raw_data_length);
HRF=HRF(1:HRF_length);

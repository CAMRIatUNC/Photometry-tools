function [HRF]=HRF(neuronal_activity,hemodynamics,sampling_rate,HRF_length,pre)
%%
% calculate HRF via simultaneously recorded neuronal_activity and
% hemodynamics time-courses.
%
% neuronal_activity: neuronal activity time-course in a single column.
% hemodynamics: hemodynamics changes time-course in a single column.
% pre: use 0 in most of the case,
%      unless you assume the recorded neuronal activity and hemodynamic 
%      response are drived by a common sourse instead of under a causal 
%      relationship from neuronal activity to hemodynamic response.
%
% Example: HRF(GCaMP_timecourse,Rhodamine_timecourse,10,25,0)
% Tzu-Hao Harry Chao 2021/10/22
%%

hemodynamics=hemodynamics-min(hemodynamics);

HRF_length=HRF_length*sampling_rate+1; % sec to data points
pre=pre*sampling_rate; % sec to data points

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

HRF=pinv(X(pre+1:raw_data_length,:))*hemodynamics(1:raw_data_length-pre);
HRF=HRF(1:HRF_length);

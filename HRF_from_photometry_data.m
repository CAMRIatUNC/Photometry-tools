function [HRF]=HRF_from_photometry_data(neuronal_activity,BOLD,sampling_rate,HRF_length,pre)
%%
% calculate HRF via simultaneously recorded neuronal_activity and BOLD time
% course.
% neuronal_activity / BOLD: data recoreded in column
% pre: 0 or any number if the neuronal_activity is considered as feedback
% response and can be after the initiation of hemodynamic response.
% Tzu-Hao Harry Chao 2020/07/04
%%

DA_zscore=(neuronal_activity-mean(neuronal_activity))/std(neuronal_activity);
BOLD=BOLD-min(BOLD);

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

HRF=pinv(X(pre+1:raw_data_length,:))*BOLD(1:raw_data_length-pre);
HRF=HRF(1:HRF_length);

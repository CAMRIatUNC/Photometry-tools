function [coef,data_fixed,Hbs,HbT]=Hb_calculation

% setpath for parameters
Rscript='/Library/Frameworks/R.framework/Versions/4.0/Resources/bin/Rscript';
Rfile='/Users/Mac/Documents/MATLAB/hemo_correction_script50/hemo_correction_script50_quad.R';
parameters_Td='~/Documents/MATLAB/hemo_correction_script50/parameters_Td1_Quad.xlsx';
load('~/Documents/MATLAB/hemo_correction_script50/References.mat','RefSpecA');

reference = RefSpecA.G_Td(:,2:3);

IntegrationTime_Sec=0.1; % Typically 10 Hz in our lab

clc
% get dataID and data path
[dataID,path_data] = uigetfile('*.txt','Select data');
cd(path_data)

% load spectra and time_stamp
class=[];
for j=1:1044
class=[class '%f '];
end
file = fopen([path_data dataID],'r');
data = textscan(file, ['%s' '%f' class],'HeaderLines',starting_row);
fclose(file);

data=cell2mat(data(2:1046));
data=(data(~isnan(data)));
data = reshape(data,length(data)./1045,1045);
timestamp_sec=data(:,1);
timestamp_sec=(timestamp_sec-timestamp_sec(1))./10^3;
data=data(:,2:1045);

disp('Fixing lost frames...')
data_fixed=zeros(length(0:IntegrationTime_Sec:timestamp_sec(end)),1044);
for i=1:1044
    data_fixed(:,i)=interp1(timestamp_sec, data(:,i),...
    0:IntegrationTime_Sec:timestamp_sec(end));
end

disp('Unmixing GCaMP and Tdtomato...')
data_fixed=data_fixed';
coef=zeros(size(reference,2),size(data_fixed,2));
for i=1:size(data_fixed,2)
coef(:,i)=max(0,lsqnonneg(reference(70:550,:),data_fixed(70:550,i)));
end
data_fixed=data_fixed';

disp('Saving xlsx for [Hb] calculation...')
COL = 290:339;
data_4Hb=data_fixed(:,COL)-(coef(1,:)'*reference(COL,1)');
t=array2table(data_4Hb);
writetable(t,[dataID(1:end-4),'.xlsx'],'WriteVariableNames' ,0);

eval(['!',Rscript,' ',Rfile,' ',dataID(1:end-4),'.xlsx ',...
    parameters_Td,' ',dataID(1:end-4),'_Hb.txt ','.'])

file = fopen([dataID(1:end-4),'_Hb.txt'],'r');
Hbs = textscan(file, ['%f' '%f' '%f' '%f'],'HeaderLines',1);
fclose(file);
Hbs=cell2mat(Hbs);
HbT=Hbs(:,1)+Hbs(:,2);

disp('Saving .mat file for all results')
save([dataID(1:end-4),'.mat'],'coef','data_fixed','Hbs','HbT')
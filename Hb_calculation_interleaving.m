function Hb_calculation_interleaving(filename,BG_filename,sampling_points,BG_sampling_points,freq)

% filename: *.txt data file exported from Oceanview
% BG_filename: background recorded w/o connecting the fiber cable to animal
% sampling_points: how many data time points to expect if there's no frame-lost
% BG_sampling_points: how many data time points to be used from the background file
% freq: simpling frequency

%%
Rscript='/Library/Frameworks/R.framework/Versions/4.0/Resources/bin/Rscript';
Rfile='/Users/Mac/Documents/MATLAB/hemo_correction_script_interleaving/hemo_correction_script400nm.R';
parameters_Td='~/Documents/MATLAB/hemo_correction_script_interleaving/parameters400.xlsx';

%%
disp('resolving interleaving spectra...')

D=dir(filename);
for data_num=1:length(D)
path_data = [D(data_num).folder,'/'];
dataID = D(data_num).name;
disp(dataID)
variableID = D(data_num).name;
variableID = variableID(1:3);
eval(['[',variableID,'_488,',variableID,'_400]','=interleaving_fix(path_data,dataID,sampling_points,freq);'])
close all
end

BG=dir(BG_filename);
path_BG = [BG(1).folder,'/'];
BG_ID = BG(1).name;
disp(BG_ID)
[BG_488,BG_400]=interleaving_fix(path_BG,BG_ID,BG_sampling_points,freq);

mean_BG_400=nanmean(BG_400,2);
mean_BG_488=nanmean(BG_488(:,2:end),2);

%%
disp('subtracting background for 400 spectra...')

grp_400=who([variableID,'*400']);
mean2D_BG_400=[];
for i=1:sampling_points
mean2D_BG_400=[mean2D_BG_400,mean_BG_400];
end
for i=1:length(grp_400)
    eval([cell2mat(grp_400(i)),'_BGfix=',cell2mat(grp_400(i)),'-mean2D_BG_400;'])
end

%%
disp('subtracting background for 488 spectra...')

grp_488=who([variableID,'*488']);

mean2D_BG_488=[];
for i=1:sampling_points
mean2D_BG_488=[mean2D_BG_488,mean_BG_488];
end
for i=1:length(grp_488)
    eval([cell2mat(grp_488(i)),'_BGfix=',cell2mat(grp_488(i)),'-mean2D_BG_488;'])
end

%%
disp('extracting 400nm time courses...')
clear grp_400 i t_400
grp_400=who('*400_BGfix');

for i=1:length(grp_400)
    eval(['ref=nanmean(',cell2mat(grp_400(i)),',2);'])
    eval(['coef=zeros(size(ref,2),size(',cell2mat(grp_400(i)),',2));'])
    for j=1:sampling_points
        eval(['coef(j)=max(0,lsqnonneg(ref(200:230,:),',cell2mat(grp_400(i)),'(200:230,j)));'])
        t_400(i,j)=coef(j);
    end
end
t_400(1)=t_400(2);



%%
disp('extracting 488nm time courses...')
clear grp_488 i t_488
grp_488=who('*488_BGfix');
for i=1:length(grp_488)
    eval(['ref=nanmean(',cell2mat(grp_488(i)),',2);'])
    eval(['coef=zeros(size(ref,2),size(',cell2mat(grp_488(i)),',2));'])
    for j=1:sampling_points
        eval(['coef(j)=max(0,lsqnonneg(ref(200:230,:),',cell2mat(grp_488(i)),'(200:230,j)));'])
        t_488(i,j)=coef(j);
    end
end
t_488(1)=t_488(2);

%% for Hb calculation
disp('calculating Hbs...')

COL = 88:152;
res_400_BGfix(:,1)=res_400_BGfix(:,2);
data_fixed=res_400_BGfix';
data_4Hb=data_fixed(:,COL);
t=array2table(data_4Hb-min(data_4Hb(:)));
writetable(t,[dataID(1:end-4),'.xlsx'],'WriteVariableNames' ,0);

eval(['!',Rscript,' ',Rfile,' ',dataID(1:end-4),'.xlsx ',...
parameters_Td,' ',dataID(1:end-4),'_Hb.txt ','.'])

file = fopen([dataID(1:end-4),'_Hb.txt'],'r');
Hbs = textscan(file, ['%f' '%f' '%f' '%f'],'HeaderLines',1);
fclose(file);
Hbs=cell2mat(Hbs);
HbT=Hbs(:,1)+Hbs(:,2);

%%
disp('saving result...')

clear BG BG_ID BG_sampling_points
clearvars -except BG* mean_BG* t_* sampling_points *_BGfix Hbs HbT dataID
save([dataID(1:end-4),'.mat'])
disp('done!!')
end

function [spec1_fixed,spec2_fixed]=interleaving_fix(path_data,dataID,sampling_points,freq)
i=0; test = {{'a'}}; % skipping headers
while isnan(str2double(test{1,1}))==1
file = fopen([path_data dataID],'r');
test = textscan(file, '%s',1,'HeaderLines',i);
fclose(file);
i=i+1;
end

class=[];
for j=1:1044
class=[class '%f '];
end
file = fopen([path_data dataID],'r');
data = textscan(file, ['%s' '%s' class],'HeaderLines',i);

timestamp=data{1,1};
timestamp_sec=[];
for i=1:length(timestamp)
test=char(timestamp(i));
if length(test)>1
timestamp_sec=[timestamp_sec;str2double(test(1:2))*3600+str2double(test(4:5))*60+str2double(test(7:end))];
end
end
timestamp_sec=timestamp_sec-timestamp_sec(1);

data=cell2mat(data(3:1046));
data=(data(~isnan(data)));
data = reshape(data,length(data)./1044,1044)';
fclose(file);
%%

AUC=sum(data(80:170,:),1); % 180:260
threshold=mean(AUC);

spec1=data(:,find(AUC<threshold));
spec1_index=find(AUC<threshold);
spec1_timestamp=timestamp_sec(find(AUC<threshold));
[spec1_timestamp, index1] = unique(spec1_timestamp);

spec2=data(:,find(AUC>threshold));
spec2_index=find(AUC>threshold);
spec2_timestamp=timestamp_sec(find(AUC>threshold));
[spec2_timestamp, index2] = unique(spec2_timestamp);

for i=1:size(spec1,1)
spec1_fixed(i,:)=interp1(spec1_timestamp,spec1(i,index1),0:1/freq:(sampling_points/freq-1/freq));
end
for i=1:size(spec2,1)
spec2_fixed(i,:)=interp1(spec2_timestamp,spec2(i,index2),0:1/freq:(sampling_points/freq-1/freq));
end

end
function interleaving_fix_withBG(filename,BG_name,dummy,formal_data_points,BG_sampling_points,freq,method,GroupOrIndividual)
%%
% 1. Read 400/488 nm interleaved recorded GCaMP spectra time series output from Oceanview.
% 2. Resolve 400 and 488 excited GCaMP spectra, and fix lost-frames if there's any.
% 3. Remove the autofluorescence of optical cable.
% 4. Calculate the amplitude time courses of GCaMP signal excited by 400 and 488 nm lasers respectively.
% 5. Extract hemoglobin changes time course from GCaMP signal excited by 400 nm laser.

% filename: filename of Oceanview recorded .txt file
% BG_name: background recored w/o connecting the fiber cable to the animal
% dummy: how many time points to skip
% formal_data_points: how many data points to expect if no frame-lost
% BG_sampling_points: how many data points to use from the background file
% freq: simpling frequency
% method: "time" or "index", fix frame-lost according to timestamps or interleaving order
% GroupOrIndividual: input "group" or "individual" (Hbs only be calculated when individual is used)

% Tzu-Hao Harry Chao 2020/06/13
%%
sampling_points=formal_data_points+dummy;

%%
disp('resolving interleaving spectra...')
D=dir(filename);
for data_num=1:length(D)
path_data = [D(data_num).folder,'/'];
dataID = D(data_num).name;
disp(dataID)
variableID = D(data_num).name;
variableID = variableID(1:3);
eval(['[',variableID,'_488,',variableID,'_400]','=interleaving_fix(path_data,dataID,sampling_points,freq,method);'])
close all
end

BG=dir(BG_name);
path_BG = [BG(1).folder,'/'];
BG_ID = BG(1).name;
disp(BG_ID)
[BG_488,BG_400]=interleaving_fix(path_BG,BG_ID,BG_sampling_points,freq,method);

mean_BG_400=nanmean(BG_400,2);
mean_BG_488=nanmean(BG_488(:,2:end),2);

%%
disp('subtracting background for 400 spectra...')
if strcmp(GroupOrIndividual,'individual')
    grp_400=who([variableID,'*400']);
elseif strcmp(GroupOrIndividual,'group')
    grp_400=who('E*400');
end
mean2D_BG_400=[];
for i=1:sampling_points
mean2D_BG_400=[mean2D_BG_400,mean_BG_400];
end
for i=1:length(grp_400)
    eval([cell2mat(grp_400(i)),'_BGfix=',cell2mat(grp_400(i)),'-mean2D_BG_400;'])
end
% for i=1:length(grp_400)
% for j=1:sampling_points
% clc
% disp(cell2mat(grp_400(i)))
% disp(j)
% eval([cell2mat(grp_400(i)),'_BGfix(:,j)=',cell2mat(grp_400(i)),'(:,j)-mean_BG_400;'])
% end
% end

%%
disp('subtracting background for 488 spectra...')
if strcmp(GroupOrIndividual,'individual')
    grp_488=who([variableID,'*488']);
elseif strcmp(GroupOrIndividual,'group')
    grp_488=who('E*488');
end
mean2D_BG_488=[];
for i=1:sampling_points
mean2D_BG_488=[mean2D_BG_488,mean_BG_488];
end
for i=1:length(grp_488)
    eval([cell2mat(grp_488(i)),'_BGfix=',cell2mat(grp_488(i)),'-mean2D_BG_488;'])
end
% for i=1:length(grp_488)
% for j=1:sampling_points
% clc
% disp(cell2mat(grp_488(i)))
% disp(j)
% eval([cell2mat(grp_488(i)),'_BGfix(:,j)=',cell2mat(grp_488(i)),'(:,j)-mean_BG_488;'])
% end
% end


%%
disp('extracting 400nm time courses...')
clear grp_400 i t_400
grp_400=who('*400_BGfix');
% for i=1:length(grp_400)
% eval(['t_400(i,:)=mean(',cell2mat(grp_400(i)),'(200:221,:));'])
% end
for i=1:length(grp_400)
    eval(['ref=nanmean(',cell2mat(grp_400(i)),',2);'])
    eval(['coef=zeros(size(ref,2),size(',cell2mat(grp_400(i)),',2));'])
    for j=1:sampling_points
        %eval(['coef(j)=max(0,lsqnonneg(ref(find(ref>(max(ref)/2)),:),',cell2mat(grp_400(i)),'(find(ref>(max(ref)/2)),j)));'])
        eval(['coef(j)=max(0,lsqnonneg(ref(200:230,:),',cell2mat(grp_400(i)),'(200:230,j)));'])
        t_400(i,j)=coef(j);
    end
end
if strcmp(GroupOrIndividual,'individual')
    t_400(1)=t_400(2);
elseif strcmp(GroupOrIndividual,'group')
    t_400(:,1)=t_400(:,2);
end


%%
disp('extracting 488nm time courses...')
clear grp_488 i t_488
grp_488=who('*488_BGfix');
% for i=1:length(grp_488)
% eval(['t_488(i,:)=mean(',cell2mat(grp_488(i)),'(200:221,:));'])
% end
for i=1:length(grp_488)
    eval(['ref=nanmean(',cell2mat(grp_488(i)),',2);'])
    eval(['coef=zeros(size(ref,2),size(',cell2mat(grp_488(i)),',2));'])
    for j=1:sampling_points
        %eval(['coef(j)=max(0,lsqnonneg(ref(find(ref>(max(ref)/2)),:),',cell2mat(grp_488(i)),'(find(ref>(max(ref)/2)),j)));'])
        eval(['coef(j)=max(0,lsqnonneg(ref(200:230,:),',cell2mat(grp_488(i)),'(200:230,j)));'])
        t_488(i,j)=coef(j);
    end
end
if strcmp(GroupOrIndividual,'individual')
    t_488(1)=t_488(2);
elseif strcmp(GroupOrIndividual,'group')
    t_488(:,1)=t_488(:,2);
end

%%
subplot(2,1,1)
if size(t_488,1)==1
plot([1:sampling_points-dummy]./10,(t_488(:,dummy+1:sampling_points)./mean(t_488(:))*100)-100,'color',[0 0.6 0])
else
plot([1:sampling_points-dummy]./10,(mean(t_488(:,dummy+1:sampling_points))./mean(t_488(:))*100)-100,'color',[0 0.6 0])
end
title('GCaMP excited by 488nm')
subplot(2,1,2)
if size(t_400,1)==1
plot([1:sampling_points-dummy]./10,(t_400(:,dummy+1:sampling_points)./mean(t_400(:))*100)-100,'color',[0.6 0 0.6])
else
plot([1:sampling_points-dummy]./10,(mean(t_400(:,dummy+1:sampling_points))./mean(t_400(:))*100)-100,'color',[0.6 0 0.6])
end
title('GCaMP excited by 400nm')

%%
if strcmp(GroupOrIndividual,'individual')
disp('calculating Hbs (applicable for spec A only)...')
Rscript='/Library/Frameworks/R.framework/Versions/4.0/Resources/bin/Rscript';
Rfile='/Users/Mac/Documents/MATLAB/hemo_correction_script_interleaving/hemo_correction_script400nm.R';
parameters_Td='~/Documents/MATLAB/hemo_correction_script_interleaving/parameters400.xlsx';
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
end
%%

disp('saving result...')
clear BG BG_ID BG_sampling_points
clearvars -except E* BG* mean_BG* t_* dummy sampling_points *_BGfix Hbs HbT dataID
save([dataID(1:end-4),'.mat'])
disp('done!!')
end



function [spec1_fixed,spec2_fixed]=interleaving_fix(path_data,dataID,sampling_points,freq,method)

% [dataID,path_data] = uigetfile('*.txt','Select data');
% cd(path_data) 
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
%figure;plot(spec1,'DisplayName','x')

spec2=data(:,find(AUC>threshold));
spec2_index=find(AUC>threshold);
spec2_timestamp=timestamp_sec(find(AUC>threshold));
[spec2_timestamp, index2] = unique(spec2_timestamp);
%figure;plot(spec2,'DisplayName','x')

if strcmp(method,'index')
    x=0;
    for i=1:length(spec1_index)-1
    space=spec1_index(i+1)-spec1_index(i);
    if space>2
    x=x+round((space+1)/2);
    else
    x=x+1;
    end
    spec1_index_corr(i)=x;
    end
    spec1_index_corr=[spec1_index_corr,x+1];

    for i=1:size(spec1,1)
    spec1_fixed(i,:)=interp1(spec1_index_corr,spec1(i,:),1:sampling_points);
    end
    
    x=0;
    for i=1:length(spec2_index)-1
    space=spec2_index(i+1)-spec2_index(i);
    if space>2
    x=x+round((space+1)/2);
    else
    x=x+1;
    end
    spec2_index_corr(i)=x;
    end
    spec2_index_corr=[spec2_index_corr,x+1];

    for i=1:size(spec2,1)
    spec2_fixed(i,:)=interp1(spec2_index_corr,spec2(i,:),1:sampling_points);
    end

elseif strcmp(method,'time')
    for i=1:size(spec1,1)
    spec1_fixed(i,:)=interp1(spec1_timestamp,spec1(i,index1),0:1/freq:(sampling_points/freq-1/freq));
    end
    
    for i=1:size(spec2,1)
    spec2_fixed(i,:)=interp1(spec2_timestamp,spec2(i,index2),0:1/freq:(sampling_points/freq-1/freq));
    end
end

end

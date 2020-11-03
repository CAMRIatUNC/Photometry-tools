function [] = unmixing2
%%
% Read raw mixed spectra time series from Oceanview output, and unmix the spectra time series one by one using user defined spectra reference.
% Tzu-Hao Harry Chao 2020/02/14
%%

clc

[dataID,path_data] = uigetfile('*.txt','Select data');
cd(path_data)

%[refID,path_ref] = uigetfile('*.csv','Select reference');
%ref=csvread([path_ref refID],1,1);

[refID,path_ref] = uigetfile('*.mat','Select reference');
ref=load([path_ref refID]);
ref=struct2cell(ref);
ref=cell2mat(ref(1));

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
data=cell2mat(data(3:1046));
data=(data(~isnan(data)));
data = reshape(data,length(data)./1044,1044)';
fclose(file);

coef=zeros(size(ref,2),size(data,2));
for i=1:size(data,2)
coef(:,i)=max(0,lsqnonneg(ref(70:550,:),data(70:550,i)));
%clc
%[num2str(i/size(data,2)*100) '%']
end


for i=1:size(coef,1)
subplot(size(coef,1),1,i)
plot(0.1:0.1:size(data,2)/10,coef(i,:)/mean(coef(i,:),2)*100-100)
%title('GCaMP time course','FontWeight','bold','FontSize',12)
xlabel('Time (s)','FontWeight','bold','FontSize',12)
ylabel('dF/F (%)','FontWeight','bold','FontSize',12)
end

% subplot(2,1,2)
% plot(0.1:0.1:size(data,2)/10,coef(2,:)/mean(coef(2,:),2)*100-100)
% title('Rhodamine time course','FontWeight','bold','FontSize',12)
% xlabel('Time (s)','FontWeight','bold','FontSize',12)
% ylabel('dF/F (%)','FontWeight','bold','FontSize',12)

save([dataID(1:length(dataID)-4) '_test.mat'],'coef')

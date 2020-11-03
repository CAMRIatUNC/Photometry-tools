function [data,timestamp_sec,wavelength]=get_photometry
%%
% read the raw txt file output from Oceanview into matlab matrix
%%
[dataID,path_data] = uigetfile('*.txt','Select data');
cd(path_data)
starting_row=0; test = {{'a'}}; % skipping headers
while isnan(str2double(test{1,1}))==1
file = fopen([path_data dataID],'r');
test = textscan(file, '%s',1,'HeaderLines',starting_row);
fclose(file);
starting_row=starting_row+1;
end

class=[];
for j=1:1044
class=[class '%f '];
end
file = fopen([path_data dataID],'r');
data = textscan(file, ['%s' '%s' class],'HeaderLines',starting_row);
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
data = reshape(data,length(data)./1044,1044);
fclose(file);

file = fopen([path_data dataID],'r');
wavelength=textscan(file, '%f','HeaderLines',starting_row-1);
wavelength=cell2mat(wavelength)';
fclose(file);

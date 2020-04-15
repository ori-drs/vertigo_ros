clc;
clear; warning off;
format longg
% Plot the sensor's timestamps to check they are increasing linearly and to
%compute the averaged frequency and its standard deviation

% Reads file existing in the home directory
% temp = readtable('manhattanOlson3500-wFalseLC.csv',);
% temp2_double = table2array(temp(:,2:5));

% plot all vertices - raw measurements
% plot(temp2_double(1:3500,2), temp2_double(1:3500,3), 'r');


ftoread = 'manhattanOlson3500-wFalseLC.csv';
fid = fopen(ftoread);
% fgetl(fid); %reads line but does nothing with it
% fgetl(fid);
% fgetl(fid);
M = textscan(fid, '%s%d%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
fclose (fid);

tail_vec =[];
head_vec =[];
wrong_lc_id_vec = [];
for i=1:length(M{1})
    str = char(M{1}(i));
    if strcmp(str,'EDGE_SE2_SWITCHABLE') ==1
        tail_vec =[tail_vec;M{2}(i)];
        head_vec =[head_vec;M{3}(i)];
        % the last n ids belong to wrong loop closures concatenated in the
        % end
        wrong_lc_id_vec = [wrong_lc_id_vec;M{2}(i)];
    end
end

plot(M{3}(1:3500),M{4}(1:3500),'r')
hold on

for i=1:length(tail_vec)
    plot([M{3}(tail_vec(i)+1) M{3}(head_vec(i)+1)], [M{4}(tail_vec(i)+1) M{4}(head_vec(i)+1)],'b')
end
% Distinct wrong loop closures
for i=0:99
    plot([M{3}(tail_vec(end-i)+1) M{3}(head_vec(end-i)+1)], [M{4}(tail_vec(end-i)+1) M{4}(head_vec(end-i)+1)],'Color', [170/255 170/255 170/255])
end









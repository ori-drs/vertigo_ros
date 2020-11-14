clc;
clear; warning off;
format longg
% Plot the slam traj and the indicate the right and wrong loop closures
% with different colours

% Reads file existing in the home directory
% temp = readtable('manhattanOlson3500-wFalseLC.csv',);
% temp2_double = table2array(temp(:,2:5));

% plot all vertices - raw measurements
% plot(temp2_double(1:3500,2), temp2_double(1:3500,3), 'r');

% insert the number of outliers, note: outliers always added at the bottom
% of the file
outliers = input('How many outliers included in the graph? '); % in long new college example, it is 116 

ftoread = 'newcollege_long.csv';
fid = fopen(ftoread);
M = textscan(fid, '%s%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); % you will need to change the number of values to match your file %f for numbers and %s for strings.
fclose (fid);

tail_vec =[];
head_vec =[];
wrong_lc_id_vec = [];
count_poses = 0;
for i=1:length(M{1})
    str = char(M{1}(i));
    if (strcmp(str,'EDGE_SE3:QUAT') ==1 && M{2}(i) ~= M{3}(i)-1)
        tail_vec =[tail_vec;M{2}(i)];
        head_vec =[head_vec;M{3}(i)];
        % the last n ids belong to wrong loop closures concatenated in the
        % end
        wrong_lc_id_vec = [wrong_lc_id_vec;M{2}(i)];
    end
    % count the number of poses in the graph
    if (strcmp(str,'VERTEX_SE3:QUAT') ==1)
       count_poses = count_poses + 1;
    end
    
end

plot(M{3}(1:count_poses),M{4}(1:count_poses),'g')
xlabel ('x (m)')
ylabel ('y (m)')
hold on
axis equal
axis tight

legend('Trajectory', 'Location', 'SouthWest')

% Plot true-positive loop closures
for i=1:length(tail_vec)-outliers
    if i==1
         plot([M{3}(tail_vec(i)+1) M{3}(head_vec(i)+1)], [M{4}(tail_vec(i)+1) M{4}(head_vec(i)+1)],'b', 'LineWidth', 3, 'DisplayName', 'True-Positive')
    end            
    h1 = plot([M{3}(tail_vec(i)+1) M{3}(head_vec(i)+1)], [M{4}(tail_vec(i)+1) M{4}(head_vec(i)+1)],'b', 'LineWidth', 3);
    set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

% Distinct wrong loop closures
for i=0:outliers-1
    if i==1
        plot([M{3}(tail_vec(end-i)+1) M{3}(head_vec(end-i)+1)], [M{4}(tail_vec(end-i)+1) M{4}(head_vec(end-i)+1)],'Color', [255/255 0/255 0/255], 'LineWidth', 3,  'DisplayName', 'False-Positive')
    end
    h2 = plot([M{3}(tail_vec(end-i)+1) M{3}(head_vec(end-i)+1)], [M{4}(tail_vec(end-i)+1) M{4}(head_vec(end-i)+1)],'Color', [255/255 0/255 0/255], 'LineWidth', 3);
    set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end







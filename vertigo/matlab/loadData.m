% This function loads CSV data from file into a timeseries object.
% INPUT:
% filename: csv input file. Format: s, ns, x, y, z, qx, qy, qz, qw.
% name    : name of data, e.g. TSIF
% start   : optional argument to specify the start of the data (as a
% fraction).
% end     : optional argument to specify the end of the data (as a
% fraction).

function ts = loadData(filename, name, start, finish)
    if (nargin<3)
      start  = 0;
      finish = 1;
    end

    raw = readtable(filename,'Delimiter',' ');    
    
    t      = raw.Var1;
    x(:,1) = raw.Var2;
    x(:,2) = raw.Var3;
    x(:,3) = raw.Var4;
    % Change quaternion to Matlab w,x,y,z format.
    q(:,1) = raw.Var8;
    q(:,2) = raw.Var5;
    q(:,3) = raw.Var6;
    q(:,4) = raw.Var7;
    
    % Filter the length of data
    a = round(length(t)*start);
    a = max(a, 1); % avoid values less than 1.
    b = round(length(t)*finish);
    t = t(a:b);
    x = x(a:b,:);
    q = q(a:b,:);
    
    % Remove any invalid values
    [row, ~] = find(isnan(t));
    t(row, :) = [];
    x(row, :) = [];
    q(row, :) = [];
    [row, ~] = find(t == 0);
    t(row, :) = [];
    x(row, :) = [];
    q(row, :) = [];
    
    ts = timeseries([x, q], t, 'name', name);
    
    fprintf('Finishing reading %d rows from %s\n', length(t), filename);
end


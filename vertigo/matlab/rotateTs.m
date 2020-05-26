% Function to rotate a trajectory around its initial point.
% INPUT:
% ts : time series of the data to rotate.
% yaw: yaw angle of rotation
% OUTPUT:
% ts_rot: rotated time series object.
% Note: this only currently rotates the position!

function ts_rot = rotateTs(ts, yaw)

ts_rot = ts;

R = [cos(yaw), -sin(yaw), 0;
     sin(yaw), cos(yaw) , 0;
     0       , 0        , 1];

% First set the origin as (0,0,0)
ts_rot.Data(:,1:3) = ts_rot.Data(:,1:3) - ts_rot.Data(1,1:3);
for i=1:length(ts_rot.Time)
   ts_rot.Data(i,1:3) = (R*(ts_rot.Data(i,1:3))')';
end

end
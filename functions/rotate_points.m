function [v_rot,r] = rotate_points(v,r)
R = [[cos(r), -sin(r)];[sin(r), cos(r)]];
v_rot = zeros(size(v));
for el = 1:size(v,2)
v_rot(:,el) = R*v(:,el);
end

end
function [bs_rot, ue_rot, ris_rot] = rotate_geometry(rot_angle, bs, ue, ris_el, plot_flag)
bs_rot = rotate_points(bs,rot_angle);
ris_rot = rotate_points(ris_el,rot_angle);
ue_rot = rotate_points(ue,rot_angle);
if plot_flag
    plot_geometry(ris_rot,ue_rot,bs_rot)
end

end
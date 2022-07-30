function my_polarplot_v(v, varargin)

Marker = {'-k'; '--k'};

if numel(varargin)
    flag_2D = varargin{1};
else
    flag_2D = false();
end

N=size(v,1);
theta=-pi:pi/1000:pi;

if ~flag_2D
    figure
    for i=1:size(v,2)
        %to_plot = -Inf*ones(1,numel(theta));
        s_v=abs(v(:,i)'*steer(N,theta,0.5)).^2;
        %to_plot = max(to_plot,s_v);
        polarplot(theta,s_v); hold on
    end
    thetalim([0, 180])
else
    figure
    subplot(1,2,1)
    
    for j=1:size(v,3)
    to_plot = -Inf*ones(1,numel(theta));
        for i=1:size(v,2)
            s_v=abs(v(:,i,j)'*steer(N,theta,0.5)).^2;
            to_plot = max(to_plot,s_v);
        end
        polarplot(theta,to_plot,Marker{j}); hold on
    end
%    rlim([0 N^2])
    subplot(1,2,2)
    for j=1:size(v,3)
    to_plot = -Inf*ones(1,numel(theta));
        for i=1:size(v,2)
            s_v=abs(conj(v(i,:,j))*steer(N,theta,0.5)).^2;
            to_plot = max(to_plot,s_v);
        end
        polarplot(theta,to_plot,Marker{j}); hold on
    end
end

%rlim([0 N^2])
% thetalim([-180 0])
% ax1=gca;
% ax1.ThetaAxisUnits='radians';
end
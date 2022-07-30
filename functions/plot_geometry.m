function [] = plot_geometry(ris_el,ue,bs, varargin)

if numel(varargin) > 0
    annotate_ue_flag = varargin{1};
end

figure
hold on
plot(ris_el(1,:),ris_el(2,:),'k*')
plot(ue(1,:),ue(2,:),'r*')
plot(bs(1,:),bs(2,:),'b*')
leg = legend('ris', 'ue', 'bs','Location','SouthEast');
x_lab = xlabel('x coord');
y_lab = ylabel('y coord');

leg.Interpreter = 'Latex';
x_lab.Interpreter = 'Latex';
y_lab.Interpreter = 'Latex';

if annotate_ue_flag
    delta = 0.5;
    for u = 1:size(ue,2)
        text(ue(1,u)+delta, ue(2,u)+delta, num2str(u))
    end
end

hold off

end
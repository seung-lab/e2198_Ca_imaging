function fig_petal(cell_info, ca_dsos, tuning_ordered, varargin)
% onoff: 0 = total, 1 = on, 2 = off, [1 2] = separate on and off


type_names={'37r','37d','37c','37v','7o','7id','7ic','7iv','2aw','63','72n'};
idx_panel=[1 1 1 1 2 2 2 2  3 4 5];
%idx_color=[1 6 2 3 1 6 2 3  3 3 3];
idx_color=[4 3 2 5 4 3 2 5  8 8 8];
type_onoff={0 0 0 0  0 0 0 0  0 1 0};
%r_lim = []
titles = {'', '', '', '37', '', '', '', '7o / 7i','2aw','63 (Ca = On only)','72n'};
fontsize = 12;



%{
n = 11;
mid = 6;
subplot(n,n, n*[0; 1; 2; 3]*[1 1 1 1] + [1;1;1;1]*[mid-2:mid+2])
polarplot(1,1)
hold on

subplot(n,n, n*[0; 1; 2; 3]*[1 1 1 1] + [1;1;1;1]*[mid-2:mid+2])
polarplot(1,1)
%}
%{
figure;

sz = 0.4;
   new_axis = axes('position',[ 0.3 0. sz sz]);
   polarplot(1,1)
   new_axis = axes('position',[ 0 0.3  sz sz]);
   polarplot(1,1)
   new_axis = axes('position',[ 0.6 0.3  sz sz]);
   polarplot(1,1)
   new_axis = axes('position',[ 0.3 0.6  sz sz]);
   polarplot(1,1)
   new_axis = axes('position',[ 0.35 0.35  0.3 0.3]);
   polarplot(1,1)
%}
idx = 1:4;
flower(type_names(idx), idx_panel(idx), idx_color(idx));
title('37')
idx = 5:8;
flower(type_names(idx), idx_panel(idx), idx_color(idx));
title('7o/7i')



function flower(type_names, idx_panel, idx_color)

figure;
sz = 0.3;
sm = 0.2;
bd = 0.05;
hx = [];
axes('position',[ 0.5-sz/2 bd sz sz]);
polarplot(0,0)
hx(4) = gca;
axes('position',[ bd 0.5-sz/2  sz sz]);
polarplot(0,0)
hx(3) = gca;
axes('position',[ 1-bd-sz 0.5-sz/2  sz sz]);
polarplot(0,0)
hx(1) = gca;
axes('position',[ 0.5-sz/2 1-bd-sz  sz sz]);
polarplot(0,0)
hx(2) = gca;
axes('position',[ 0.5-sm/2 0.5-sm/2  sm sm]);
polarplot(0,0)
hx(5) = gca;

tmp = gca();
colors = tmp.ColorOrder;
colors(end+1,:) = [0 0 0];

petalColors = distinguishable_colors(30);

layer = 3;
for j=1:numel(type_names)

	cells = get_cell_info(cell_info, type_names{j});
    cell_ids = [cells.cell_id];

	for i=1:numel(cell_ids)

        idx=find(cell_ids(i)==ca_dsos.omni_id);
        if isempty(idx)
            continue;
        else
        	ca_id = ca_dsos.ca_id(idx);
        end

        % petal
        axes(hx(j));
	    theta = pi/4 * [1:8].';
	    theta = pi/2+theta;  % to final "standard" coord
	    rho = tuning_ordered(layer,:,ca_id);
	    rho = sum(rho,1);
	    polarplot(theta([end 1:end]), rho([end 1:end]),'LineWidth',1,'Color',petalColors(i,:));
        hold on;

        % middle plot
		axes(hx(5));

        theta = ca_dsos.ds_theta(idx,layer);
        rho = ca_dsos.ds_r(idx,layer);
        
        theta=pi/2+theta;  % to final "standard" coord
        polarplot([0 theta],[0 rho],'LineWidth',1,'Color',colors(idx_color(j),:));    
        hold on;
    end
end

for j=1:5
	axes(hx(j));
	applystyle(gca);
end

axes(hx(5));
ax = gca();
lim = rlim();
if lim(2)>6
    rlim([0 15])
    ax.RTick = [5 10 15];
else
    rlim([0 6])
    ax.RTick = [6];
end

figure_size_x2();
    h = gcf();
    h.Position(3) = h.Position(4);

end %nested func
end %func




function applystyle(ax)
	fontsize = 12;

	ax.FontSize = fontsize;
    ax.ThetaTick=0:45:315;
    ax.GridAlpha = 1;
    ax.GridColor = 0.8*[1 1 1];
    ax.ThetaTickLabel = {};
    ax.GridLineStyle = '--';

    lim = rlim();
        if lim(2)>=30
            rlim([0 60])
            ax.RTick(1) = []; % remove 0 tick
        else
            rlim([0 20])
            ax.RTick(1) = []; % remove 0 tick
        end

end

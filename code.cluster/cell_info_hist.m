
function [cells,cell_stat,ctype,bin]=cell_info_hist(cell_info,type_names, stat_type, x_lim, varargin) %,  p, pminus, printfigure

% stat_type: prcntile, prcntileDiff, peakWidth

nvarargin = length(varargin);
optargs = {0.2, [], [], '', [], Inf, false, false};
optargs(1:nvarargin) = varargin;
[p, pminus, divisions, printfigure, binsize, cutoff, printcells, showstrat] = optargs{:};
if isempty(cutoff)
    cutoff = Inf;
end


binstep=0;

cells=[];
for j=1:numel(type_names)
    idx=strncmp({cell_info.type},type_names{j}, length(type_names{j}));
    if isempty(find(idx))
        error(sprintf('Unrecognized type "%s"', type_names{j}));
    end
    cells=[cells; [cell_info(idx).cell_id]'];
end    
N=numel(cells);

if N==0
    error('no cells found')
end

strat=cell_info_bin_strat(cell_info,binstep);
cell_stat=zeros(size(cells));
ctype=cell(size(cells));

for j=1:N
    cell_info_elem = get_cell_info(cell_info, cells(j));

    ctype{j}=cell_info([cell_info.cell_id]==cells(j)).type;

    s=strat{cells(j)}(:,2);
    x=strat{cells(j)}(:,1);
    %%{
    s=s(x<cutoff);
    x=x(x<cutoff);
    %}

    switch stat_type
    case {'prcntile', 'ptile', 'ptileDiff', 'percntileDiff', 'ptile-'}
        %[locm,locl,locr]=cell_info_get_cell_height_from_prcntile([x s],p);
        %cell_stat(j)=locr;  %20
        %cell_stat(j)=locl; %80
        %cell_stat(j)=locl-locr;
        %tmp = sortrows([x s], 1);   % put into ascending order
        cell_stat(j) = get_percentile([x s],p);
        if ~isempty(pminus) || strcmp(stat_type, 'percntileDiff') || strcmp(stat_type, 'ptileDiff')	...
                    || strcmp(stat_type, 'ptile-')  % compute height between p and pminus
        	%[~,~,minus]=cell_info_get_cell_height_from_prcntile([x s],pminus);
            minus = get_percentile([x s],pminus);
        	cell_stat(j)=cell_stat(j) - minus;
        	if cell_stat(j) < 0
        		cell_stat(j) = -cell_stat(j);
        	end
        end

    case {'ptile+'}
        cell_stat(j) = get_percentile([x s],p);
        if ~isempty(pminus) || strcmp(stat_type, 'percntileDiff') || strcmp(stat_type, 'ptileDiff') % compute height between p and pminus
            %[~,~,minus]=cell_info_get_cell_height_from_prcntile([x s],pminus);
            minus = get_percentile([x s],pminus);
            cell_stat(j)=cell_stat(j) + minus;
        end

    case {'peak_no_split', 'PNS'}
        [x1, x2]=cell_info_get_cell_height_from_peak_not_allow_split([x s],p);
        cell_stat(j) = x2 - x1;

    case {'peaks', 'PS'}
        [x1, x2]=cell_info_get_cell_height_from_peak_allow_split([x s],p);
        %cell_stat(j) = x2 - x1;
        xx1=min([x1; x2]);
        xx2=max([x1; x2]);
        cell_stat(j) = xx2 - xx1;

    case {'asymmetry', 'asym'}
        cell_stat(j) = norm(cell_info_elem.asymm_index);

    case {'asym_2an', 'asym2'}
        cell_stat(j) = log(norm(cell_info_elem.asymm_2an));
        if (cell_info_elem.is_cutoff)
            cell_stat(j) = NaN;
        end
    case {'asym_2an_p', 'asym2_p'}
        cell_stat(j) = log(norm(cell_info_elem.asymm_2an_prj));
        if (cell_info_elem.is_cutoff)
            cell_stat(j) = NaN;
        end

    case {'density'}
        cell_stat(j) = cell_info_elem.area_projection / cell_info_elem.area_hull;

    case {'size'}
        cell_stat(j) = cell_info_elem.area_hull;

    case {'maxdiameter'}
        cell_stat(j) = cell_info_elem.max_diameter;

    case {'radius'}
        cell_stat(j) = cell_info_elem.radii_hull;

    case {'area'}
        cell_stat(j) = cell_info_elem.area_projection;

    case {'value'}
        cell_stat(j) = mean( s(find(x>p & x<p+1)) );

    case {'slope'}
        xx = x(end:-1:1);
        yy = s(end:-1:1);
        %xx = xx(5<x<30);
        %yy = yy(5<x<30);
        %xx = x(5<x<25);
        %yy = s(5<x<25);
        slope = [xx(:) ones(length(xx),1)] \ yy(:);
        cell_stat(j) = slope(2);
        %yy = yy / sum(yy);
        yy = yy - mean(yy);
        half = (length(yy)-1)/2;
        slope = [-half:half] * yy(:);
        slope = yy(find(xx>25, 1, 'first')) - yy(find(xx>5, 1, 'first'));
        cell_stat(j) = slope;

    case {'entropy'}
        % test s = cell_info_elem.strat_unrml(:,2);
        % Renormalize to make sure?
        e = log(s);
        e(s==0) = 0;
        e = - s(:).' * e(:);
        cell_stat(j) = e;

    case {'s:t.stratratio'}
        as = sum(s(x<28)) + sum(s(x>62));
        at = sum(s(x>28 & x<62));
        if at
            cell_stat(j) = log(as/at);
        else
            cell_stat(j) = 20;
        end

    case {'on:off.stratratio'}
        on = sum(s(x>45));
        off = sum(s(x<45));
        if on && off
            cell_stat(j) = log(on/off);
        elseif off
            cell_stat(j) = -20 
        else
            cell_stat(j) = 20;
        end
%{
    case {'on-off'}
        on = sum(s(x>45));
        off = sum(s(x<45));
        cell_stat(j) = on-off;

    case {'trans_on-trans_off'}
        cell_stat(j) = cell_info_get_strat_property(cell_info_elem, 'trans_on-trans_off');

    case {'sus-trans+sus_on-trans_on'}
        cell_stat(j) = cell_info_get_strat_property(cell_info_elem, 'sus-trans') ...
                     + cell_info_get_strat_property(cell_info_elem, 'sus_on-trans_on');
%}
    otherwise
        %try
            cell_stat(j) = cell_info_get_strat_property(cell_info_elem, stat_type);
        %catch
        %error('not recognized stat name') 
        %end

    end
end
%cell_stat

autoxlim = isempty(x_lim);
if ~autoxlim && ( x_lim(2) > 20 * max(cell_stat) || x_lim(2) < min(cell_stat) )
    autoxlim = true;
    warning('Auto resetting improper limits');
end
if autoxlim
    x_lim(1) = min(cell_stat);
    x_lim(2) = max(cell_stat);
end

range = x_lim(2)-x_lim(1);
if ~isempty(binsize) && ( binsize >= range || binsize < range/1000 )
    binsize = [];
    warning('Auto resetting improper binsize');
end
if isempty(binsize)
    if range < 100 && range > 5
        binsize = 1;
    else
        binsize = range / 20;
    end
end

ub = x_lim(2);
if mod(range, binsize) ~= 0
    % fix missing last bin % TODO: maybe I should use the native support for automatic binning in the histcounts().
    ub = ub + binsize;
end
binranges=x_lim(1):binsize:ub;
%[cnts,bin]=histc(cell_stat,binranges);
%bar(binranges,cnts,'histc');
cnts = histcounts(cell_stat,binranges);
bin = discretize(cell_stat,binranges);
%binranges
%cnts = histcounts(cell_stat(~isnan(cell_stat)), binranges);
%bin = discretize(cell_stat(~isnan(cell_stat)), binranges);
figure;
bar(binranges(1:end-1),cnts,'histc');


strat_mean = [];
%%{
% order by stat within each type
for type = type_names(:).'
    type = type{1};
    idx = find(strncmp(ctype, type, length(type)));
    [B,I] = sort(cell_stat(idx));
    cells(idx)     = cells(idx(I));
    cell_stat(idx) = cell_stat(idx(I));
    bin(idx)       = bin(idx(I));
    tmp = cat(3, strat{cells(idx)});
    strat_mean(:, end+1) = squeeze(mean(tmp(:,2,:), 3));
end
%}



ax=gca();
if ~autoxlim
    ax.XLim=x_lim;
end

%ax.XTick=binranges;
%ax.XTick=0:2:50;

%%{
if length(type_names) > 1
    ymax = ax.YLim(2);
    ystep = ymax / 15;
    y = ymax/2;
    if y < 10.6
        y = 10.6;
    end

    linesx = [];
    linesy = [];
    for type = type_names(:).'
        type = type{1};
        stats = cell_stat(strncmp(ctype,type,length(type)));
        xmin = min(stats);
        xmax = max(stats);
        text(mean([xmin,xmax]), y, type); %, 'VerticalAlignment', 'bottom');
        line([xmin; xmax], [y; y]);
        line(stats(:).', repmat(y,length(stats),1), 'Marker', '.');
        %line(stats(:).', repmat(y,length(stats),1), 'Marker', '+');
        linesx = [linesx, [xmin; xmax]];
        linesy = [linesy, [y; y]];
        y = y+1;%ystep;
    end
    %line(linesx, linesy);
end
%}

%title([stat_type, '  ', strjoin(type_names)])
title([stat_type, '  ', sprintf('%g %g', p, pminus)])

if ~isempty(divisions) && max(divisions) < x_lim(2) && min(divisions) > x_lim(1)
    n = length(divisions);
    divisions = divisions(:).'; %make row vec
    %line([divisions; divisions], repmat([0; 10], 1, n), 'Color', [0 .6 .7])
    line([divisions; divisions], repmat([0; ax.YLim(2)], 1, n), 'Color', [1 .3 .7])
    xtick = ax.XTick;
    ax.XTick = unique([xtick divisions]);
end

if showstrat
    new_axis = axes('position',[ 0.6 0.7 0.4 0.3]);
    x=strat{cells(1)}(:,1);
    plot(x, strat_mean);
    set(new_axis,'color','none');
    set(new_axis,'visible','off');
    legend(type_names(:));
end


if ~isempty(printfigure) && ischar(printfigure)
    print(gcf, '-r300', printfigure, '-dpng');
    %print(gcf, '-r300', printfigure, '-depsc');
    fprintf('saved file %s \n', printfigure)
    %close(summary_fig_h);
end

if printcells
    tab = repmat(sprintf('\t'), N, 1);
    [num2str(cells) tab num2str(cell_stat) tab char(ctype)]
    %{
    for j=1:N
        fprintf(stat)
    end
    %}
end

end
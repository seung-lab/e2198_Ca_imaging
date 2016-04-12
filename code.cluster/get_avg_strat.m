function [s] = get_avg_strat(cell_info, typename)

stratname = 'strat_nrml';

sac = get_cell_info(cell_info, typename);
strat = cat(3, sac.(stratname));
strat = mean(strat, 3);
s=strat(:,2);
x=strat(:,1);
onsac = s;

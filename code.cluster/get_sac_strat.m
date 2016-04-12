function [onsac, offsac] = get_sac_strat(cell_info)

stratname = 'strat_nrml';

sac = get_cell_info(cell_info, 'ON SAC');
strat = cat(3, sac.(stratname));
strat = mean(strat, 3);
s=strat(:,2);
x=strat(:,1);
s(x>85) = 0;
onsac = s;

sac = get_cell_info(cell_info, 'OFF SAC');
strat = cat(3, sac.(stratname));
strat = mean(strat, 3);
s=strat(:,2);
x=strat(:,1);
s(x<10) = 0;
offsac = s;

%10/85
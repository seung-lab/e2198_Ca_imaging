%get_super_types.m
% GC's only
 


type_ds_sac = {'37' '7i' '7o'};

type_sus_off = {'1ni' '1no' '1ws' '1wt' '25' '27' '2an' '2aw' '2i' '2o' '3i' '3o'    '28'};
type_trans_off = {'4i' '4on' '4ow' '5to' };
type_trans_onoff = {'51' '5si' '5so' '5ti' '63'};
type_trans_on = {'6sn' '6sw' '6t'};
type_sus_on = {'28'    '58' '72n' '72w' '8' '9'}; %'81i' '81o' '82n' '82wi' '82wo' '83' '8n' '8w'  %'91n' '91w' '9m' '9n' '9w'};


types_alpha = {'1wt' '4ow' '6sw' '8w'};
types_mini_alpha = {'4i' '4on' '6sn' '8n'};
types_alpha_all = [types_alpha  types_mini_alpha];


types_oodsgc = {'37c' '37v' '37r' '37d'};
types_ondsgc = {   '7o'   '7iv' '7ir' '7id'};
types_ds_sac = [types_oodsgc  types_ondsgc];
types_ds_sac__grouped = {'37' '7i' '7o'};

types_sus_off__no_on = {'1ni' '1no'   '1wt' '2an' '2aw' '2i' '2o'   '3o'};
types_sus_off__has_on = {'1ws' '3i'   '25' '27' '28'};
types_sus_off = [types_sus_off__no_on  types_sus_off__has_on];

types_trans_off = {'4i' '4on' '4ow' };
types_trans_onoff = {'51' '5si' '5so' '5ti'  '5to' '63'};
types_trans_on = {'6sn' '6sw' '6t'};

types_sus_on__ca_slow = {'73' '72' '82wi' '82wo'  '9n' '9w'};
types_sus_on__ca_fast = {
    '81i' ...
    '81o' ...
    '82n' ...
    '85' ...
    '8n' ...
    '8w' ...
    '91-' ...
    '915' ...
};
types_sus_on = [types_sus_on__ca_slow  types_sus_on__ca_fast];

types_ca_only_off = {'1ni' '1no' '1wt' '2an' '2aw' '2i' '2o'   '3o' '4i' '4on' '4ow'};
types_ca_has_on = {'1ws' '3i' ...
    '25' ...
    '27' ...
    '28' ...
    '37c' ...
    '37d' ...
    '37r' ...
    '37v' ...
    '51' ...
    '85' ...
    '5si' ...
    '5so' ...
    '5ti' ...
    '5to' ...
    '63' ...
    '6sn' ...
    '6sw' ...
    '6t' ...
    '73' ...
    '72' ...
    '7ir' ...
    '7iv' ...
    '7id' ...
    '7o' ...
    '81i' ...
    '81o' ...
    '82n' ...
    '82wi' ...
    '82wo' ...
    '915' ...
    '8n' ...
    '8w' ...
    '9n' '9w' ...
    '91-' ...
};

assert(isequal(sort(types_ca_has_on), sort([types_ds_sac  types_sus_off__has_on  types_trans_onoff  types_trans_on  types_sus_on])))
assert(isequal(sort(types_ca_only_off), sort([types_sus_off__no_on  types_trans_off])))

ca_groups = {
    'outer marginal'            types_sus_off 
    'outer central'            type_trans_off %types_trans_off
    'inner-outer central'      type_trans_onoff  %types_trans_onoff 
    'inner central'            types_trans_on 
    'inner marginal'            types_sus_on 
};
ca_groups = cell2table(ca_groups, 'VariableNames', {'name', 'types'});

ca_groups3 = {
    'outer extra-SAC'            types_sus_off      -7
    'outer intra-SAC'            type_trans_off  -7 %types_trans_off
    'inner-outer intra-SAC'      type_trans_onoff  -7 %types_trans_onoff 
    'inner intra-SAC'            types_trans_on 0
    'inner extra-SAC: slow'      types_sus_on__ca_slow 0
    'inner extra-SAC: fast'      types_sus_on__ca_fast 0
};
ca_groups3 = cell2table(ca_groups3, 'VariableNames', {'name', 'types', 'shift'});


ca_groups2 = {
	'classical DS'               types_ds_sac           1 1 'on/off'     {'p', 'r'} {'p', 'y'} {'p', 'y'} {'p', 'y'}
    'outer extra-SAC: off only'  types_sus_off__no_on   0 1 'off'        {'>', 'm'} {'^', 'g'} {'^', 'g'} {'^', 'g'}
    'outer extra-SAC: on/off'    types_sus_off__has_on  1 1 'on/off'     {'h', 'k'} {'h', 'g'} {'h', 'g'} {'h', 'g'}
    'outer intra-SAC'            type_trans_off        0 1 'off'        {'<', 'm'} {'^', 'r'} {'^', 'r'} {'^', 'r'}
    'inner-outer intra-SAC'      type_trans_onoff      1 1 'on/off'     {'d', 'g'} {'d', 'r'} {'d', 'filled', 'MarkerFaceColor', [1 .7 .7]} {'d', 'r'}
    'inner intra-SAC'            types_trans_on         1 0 'on'         {'<', 'b'} {'v', 'r'} {'filled','vr'} {'v', 'r'}
    'inner extra-SAC: slow'      types_sus_on__ca_slow  1 NaN 'on/(off)' {'c'}      {'g'}      {'filled','g'}      {'g'}     
    'inner extra-SAC: fast'      types_sus_on__ca_fast  1 NaN 'on/(off)' {'>', 'b'} {'v', 'g'} {'filled','v','g'} {'v', 'g'}

    'classical DS'               types_ds_sac           0 1 'on/off'     {'p', 'b'} {'p', 'y'} {'p', 'b'} {'p', 'b'}
    'outer extra-SAC: on/off'    types_sus_off__has_on  0 1 'on/off'     {'h', 'b'} {'h', 'g'} {'h', 'b'} {'h', 'b'}
    'inner-outer intra-SAC'      type_trans_onoff      0 1 'on/off'     {'d', 'b'} {'d', 'r'} {'d', 'filled', 'MarkerFaceColor', 'b'} {'d', 'b'}
};
ca_groups2 = cell2table(ca_groups2, 'VariableNames', {'name', 'types', 'has_on', 'has_off', 'onoff_text', 'style1', 'style2', 'style3', 'style4'});


co = get(groot, 'defaultAxesColorOrder');  % blue, red, yellow/orange, purple, green
ca_groups4 = {
    %'classical DS'               types_ds_sac           1 1 'on/off'     {'p', 'r'} {'p', 'y'} {'filled', 'filled', 'MarkerFaceColor', co(6,:)} {'p', 'y'}
    %'outer marginal'  types_sus_off   0 1 'off'        {'>', 'm'} {'^', 'g'} {'filled', 'filled', 'MarkerFaceColor', co(1,:)} {'^', 'g'}
    %%{
    'outer marginal'  types_sus_off__no_on   0 1 'off'        {'>', 'm'} {'^', 'g'} {'filled', 'filled', 'MarkerFaceColor', co(1,:)} {'^', 'g'}
    %'outer marginal'    types_sus_off__has_on  1 1 'on/off'     {'h', 'b'} {'h', 'g'} {'filled', 'MarkerFaceColor', co(1,:)} {'h', 'b'}
    'outer marginal'    types_sus_off__has_on  0 1 'on/off'     {'h', 'b'} {'h', 'g'} {'filled', 'MarkerFaceColor', co(1,:)} {'h', 'b'}
    %}
    'outer central'            type_trans_off        0 1 'off'        {'<', 'm'} {'^', 'r'} {'filled', 'MarkerFaceColor', co(2,:)} {'^', 'r'}
    'inner-outer central'      type_trans_onoff      1 1 'on/off'     {'d', 'g'} {'d', 'r'} {'filled', 'MarkerFaceColor', co(3,:)} {'d', 'r'}
    'inner central'            types_trans_on         1 0 'on'         {'<', 'b'} {'v', 'r'} {'filled', 'MarkerFaceColor', co(4,:)} {'v', 'r'}
    'inner marginal'      types_sus_on__ca_slow  1 NaN 'on/(off)' {'c'}      {'g'}      {'filled', 'MarkerFaceColor', co(5,:)}      {'g'}     
    'inner marginal'      types_sus_on__ca_fast  1 NaN 'on/(off)' {'>', 'b'} {'v', 'g'} {'filled', 'MarkerFaceColor', co(5,:)} {'v', 'g'}

    %{
    %'inner-outer central'      type_trans_onoff      0 1 'on/off'     {'d', 'b'} {'d', 'r'} {'d', 'filled', 'MarkerFaceColor', 'b'} {'d', 'b'}
    'inner-outer central'      {'63'}      0 1 'on/off'     {'d', 'b'} {'d', 'r'} {'filled', 'MarkerFaceColor', co(3,:)} {'d', 'b'}
    'classical DS'               types_ds_sac           0 1 'on/off'     {'p', 'b'} {'p', 'y'} {'filled', 'filled', 'MarkerFaceColor', co(7,:)} {'p', 'b'}
    %}
};
ca_groups4 = cell2table(ca_groups4, 'VariableNames', {'name', 'types', 'has_on', 'has_off', 'onoff_text', 'style1', 'style2', 'style3', 'style4'});

ca_groups_alphas = {
    'outer marginal'  {{'1wt'}}   0 1 'off'        {'>', 'm'} {'^', 'g'} {'filled', 'filled', 'MarkerFaceColor', co(1,:)} {'^', 'g'}
    'outer central'            {{'4ow'}}        0 1 'off'        {'<', 'm'} {'^', 'r'} {'filled', 'MarkerFaceColor', co(2,:)} {'^', 'r'}
    'inner central'            {{'6sw'}}         1 0 'on'         {'<', 'b'} {'v', 'r'} {'filled', 'MarkerFaceColor', co(4,:)} {'v', 'r'}
    'inner marginal'      {{'8w'}}  1 NaN 'on/(off)' {'c'}      {'g'}      {'filled', 'MarkerFaceColor', co(5,:)}      {'g'}
};
ca_groups_alphas = cell2table(ca_groups_alphas, 'VariableNames', {'name', 'types', 'has_on', 'has_off', 'onoff_text', 'style1', 'style2', 'style3', 'style4'});


super_types = containers.Map();
%super_types('ds_sac') = type_ds_sac;
super_types('sus_off') = type_sus_off;
super_types('trans_off') = type_trans_off;
super_types('trans_onoff') = type_trans_onoff;
super_types('trans_on') = type_trans_on;
super_types('sus_on') = type_sus_on;


% a table is probably easier to use
%super_types.values
%super_types.keys

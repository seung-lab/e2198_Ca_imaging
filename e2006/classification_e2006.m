>> kn_e2006_ALLSKELETONS_FINAL2012_cellIDs_sortedByType_MAR2013(kn_e2006_ALLSKELETONS_FINAL2012_cellIDs_sortedByTyp_MAR2013_typp==216)
Warning: 'kn_e2006_ALLSKELETONS_FINAL2012_cellIDs_sortedByTyp_MAR2013_typp' exceeds the MATLAB maximum
name length of 63 characters and will be truncated to
'kn_e2006_ALLSKELETONS_FINAL2012_cellIDs_sortedByTyp_MAR2013_typ'. 

ans =

    73    66    82    72     5    53

  64×1 cell array

    'amfoff_1324'
    'amfoff_1419'
    'amfoff_1555'
    'amfoff_1566'
    'amfoff_1610'
    'amfoff_1644'
    'amfoff_1838'
    'amfoff_1912'
    'amfoff_723'
    'amfoff_727'
    'amfoff_826'
    'amfoff_832'
    'amfoff_844'
    'amfoff_912'
    'amfoff_936'
    'amfon_138'
    'amfon_150'
    'amfon_156'
    'amfon_259'
    'amfon_269'
    'amfon_297'
    'amfon_307'
    'amfon_370'
    'amfon_43'
    'amfon_44'
    'amfon_55'
    'amfon_65'
    'amfon_69'
    'gac_105'
    'gac_115'
    'gac_131'
    'gac_175'
    'gac_194'
    'gac_199'
    'gac_234'
    'gac_277'
    'gac_299'
    'gac_338'
    'gac_35'
    'gac_354'
    'gac_361'
    'gac_403'
    'gac_410'
    'gac_422'
    'gac_436'
    'gac_441'
    'gac_457'
    'gac_461'
    'gac_478'
    'gac_505'
    'gac_508'
    'gac_512'
    'gac_518'
    'gac_83'
    'gac_93'
    'gal_471'
    'gal_499'
    'gal_524'
    'gao_27'
    'gao_313'
    'gao_390'
    'gao_424'
    'gao_433'
    'gao_488'


figure;cell_info_scatter_plot(e2006, {'e2006gc'}, {'sac_corr' 'sac_corr'})
%e2006 = update_e2006_type(e2006, 'dsgc', {''})


% a-2
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'e2006gc' },{'sus-trans'},[-1,1],[3],[5],'Split a-2','Marginal - central',[],[],[0.05],'',[],[], 1);

e2006 = update_e2006_type(e2006, 'marginal', {
'gac_105'
'gac_518'
'gac_131'
'gac_461'
'gac_115'
'gac_93'
'gac_441'
'gac_457'
'gac_299'
'gac_277'
'gac_361'
'gac_410'
'gac_508'
'gac_175'
});


% a-3
e2006tran = [];
for elem = e2006(:).'
	elem.strat_nrml(elem.strat_nrml(:,1)<28 | elem.strat_nrml(:,1)>62, 2) = 0;
	e2006tran = [e2006tran; elem];
end
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006tran,{'e2006gc'},{'on-off'},[-1,1],[1],[5],'Split a-3','Central inner - outer',[],[],[-0.5,0.5],'',[],[],[],[]);
e2006 = update_e2006_type(e2006, 'outer central', {
'gac_505'
'gac_234'
'gac_403'
'gac_512'
});

% b-1  fail, no 5to, wrong division line location
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'outer central'},{'ptile'},[0.21 0.36],[0.5],[1],'Split b-1','10th percentile',[0.1],[],[0.28],'',[0.01],[0.9],[1],[1]);

% b-2  soma size for 4ow (likely gac_234, 505)

% b-3 fail. couldn't tell, 
%[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'outer central'},{'ptile-'},[0.05 0.11],[0.5],[1],'Split b-3','70th percentile - 15th percentile',[0.7],[0.15],[0.075],'',[0.004],[],[],[1]);
%[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'outer central'},{'ptile-'},[0.05 0.11],[0.5],[1],'Split b-3','70th percentile - 15th percentile',[0.7],[0.15],[0.075],'',[0.004],[0.8],[1],[]);

% change 0.15 to >=0.25
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'outer central'},{'ptile-'},[0.0 0.15],[0.5],[1],'Split b-3','70th percentile - 35th percentile',[0.75],[0.35],[0.075],'',[0.004],[0.8],[1],[]);
draw_e2006_cells(e2006, 'outer central')

% c-1
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'e2006gc'},{'corr','OFF SAC'},[0 0.8],[1],[5],'Split c-1','Off SAC similarity',[],[],[0.42],'',[0.05],[0.9],[1],[]);
e2006 = update_e2006_type(e2006, '63', {
    %'gac_478'
    'gac_436'
    'gac_354'
});

% c-2
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'e2006gc'},{'ptile'},[-0.05 0.45],[1],[5],'Split c-2','5th percentile',[0.05],[],[0.285],'',[0.02],[],[1],[0]);
e2006 = update_e2006_type(e2006, '5s', {
    'gac_422'
    'gac_338'
    'gac_35'
    'gac_199'
    'gac_83'
    'gac_194'
});    %10044 gac_422 wierd (type 5 in Moritz)

% c-3
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'e2006gc'},{'branch'},[0.055 0.093],[1],[5],'Split c-3','Arbor complexity (\mum^{-1})',[0.15],[],[0.0745],'',[0.002],[],[]);



% a-4  (TODO: used full strat. Only-marginal currently seems worse, w/o soma cutoff tho)
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'marginal' },{'on-off'},[-1,1],[2],[5],'Split a-4','Marginal inner - outer',[],[],[-0.15],'',[],[],1,[]);

e2006 = update_e2006_type(e2006, 'inner marginal', {
'gac_461'
'gac_299'
'gac_361'
'gac_508'
'gac_175'
'gac_410'
});


e2006sus = [];
for elem = e2006(:).'
	elem.strat_nrml(elem.strat_nrml(:,1)>28 & elem.strat_nrml(:,1)<62, 2) = 0;
	e2006sus = [e2006sus; elem];
end
% e-1
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'marginal'},{'ptile'},[0 0.9],[3],[5],'Split e-1','Marginal 80 percentile',[0.80],[],[0.45],'',[],[],[1],[0]);
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006sus,{'marginal'},{'ptile'},[0 0.9],[3],[5],'Split e-1','Marginal 80 percentile',[0.80],[],[0.45],'',[],[0.9],[1],[0]);
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006sus,{'marginal'},{'ptile'},[0 0.9],[3],[5],'Split e-1','Marginal 80 percentile',[0.85],[],[0.45],'',[],[0.89],[1],[0]);
e2006 = update_e2006_type(e2006, '27/28', {
'gac_115'
});

% e-2 soma size. not done. possibly 1wt gac_131, gac_105 2o.  gac_131 definitely 1wt

% e-3
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'marginal'},{'ptile'},[0 0.3],[1],[5],'Split e-3','50th percentile',[0.5],[],[0.11],'',[],[0.9],[1],[]);
e2006 = update_e2006_type(e2006, '1ws/1ni/1no', {
'gac_277'
});
% note (not 1ws by size)


% e-4 manual, arbor complexity
draw_e2006_cells(e2006, 'marginal')
e2006 = update_e2006_type(e2006, 'marg complex', {
'gac_105'	% 2an, possibly 2o
'gac_518'
'gac_93'
});



% e-10
% cutting off soma/axon at 0.8
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'marg complex'},{'ptile-'},[0.135 0.38],[0.5],[1],'Split e-10','90th percentile - 45th percentile',[0.9],[0.45],[0.3], '', [0.01],[0.8],[1],[0]);
e2006 = update_e2006_type(e2006, {'2an' '25' '25'}, {
'gac_105'	% 2an, possibly 2o
'gac_518'
'gac_93'
});

remaing TODO:
131 (1wt)
441
457

Type 3: 2an, 1wt


% f-1  soma size

% f-2
% cutting off soma/axon at 1 (anything >=0.9 is fine)
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'inner marginal'},{'ptile'},[0.71 0.97],[1],[5],'Split f-2','95th percentile',[0.95],[],[0.81],'',[0.015],[1],[1],[2]);
e2006 = update_e2006_type(e2006, {'bistrat inner marginal'}, {
'gac_461'   % likely 72
});


% f-6
% cutting off soma/axon 
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'inner marginal'},{'ptile'},[0.67 0.735],[1],[5],'Split f-6','80th percentile',[0.8],[],[0.71],'',[0.005],[1],[1],[2]);

% f-9
% cutting off soma/axon 
[cell_ids,stat,ctype,bin]=cell_info_hist(e2006,{'inner marginal'},{'trans'},[0.11 0.41],[1],[5],'Split f-9','Central',[],[],[0.23],'',[0.03],[1],[1],[1]);
making 'gac_461' borderline 73


'gac_361'   looks like 91/82n, strat looks like 91

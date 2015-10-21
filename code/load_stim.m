function [stim_struct] = f(full_filename)

% get handles of figs

fid=fopen(full_filename);
while 1
    tline = fgetl(fid);
    switch 1
        case ~isempty(strfind(tline, 'nConditions'))
                C = textscan(tline,'%[^=]%s%d8');
                nconds = C{3};
        case ~isempty(strfind(tline, 'CondNames'))  
                C = textscan(tline,'%[^=]%s%s');
                conds = char(C{3});
                D = textscan(conds,'%s',nconds,'delimiter',',');
                condnames = [];
                for temp = 1:length(D{1})
                    condnames = [condnames D{1}(temp)];
                end
                condnames;
        case ~isempty(strfind(tline, 'tPre'))
                C = textscan(tline,'%[^=]%s%f');
                t_pre = C{3};
        case ~isempty(strfind(tline, 'tStim'))
                C = textscan(tline,'%[^=]%s%f');
                t_stim = C{3};
        case ~isempty(strfind(tline, 'tPost'))
                C = textscan(tline,'%[^=]%s%f');
                t_post = C{3};
    end
    if ~ischar(tline),   break,   end
end
fclose(fid);

stim_struct = struct('nconds',nconds,...
                     'condnames',condnames,...
                     't_pre',t_pre,...
                     't_stim',t_stim,...
                     't_post',t_post);
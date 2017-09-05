function txt = table_row_text(row)
% print a table row to a text string

txt = strsplit(evalc('disp(row)'), '\n');
txt = txt{end-1};

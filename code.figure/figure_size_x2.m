function figure_size_x2(multiplier) 
	% Default multiplier = 2
	if ~exist('multiplier', 'var')
		multiplier = 2;
	end
    h = gcf();
    pos = h.Position;
    pos(3:4) = pos(3:4) .* multiplier;
    h.Position = pos;
end

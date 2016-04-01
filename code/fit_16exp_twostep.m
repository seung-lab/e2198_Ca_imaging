%coeffs16exp_twostep = NaN(634, 18);
%fit16exp_twostep = NaN(31*8, 634);


for ind = 90 %[346 90] %c5ti

	tauguess = 40;

	%yactual = roi_sums_means(:, :, ind);
	%yactual = yactual(:);
	yactual = roi_sums_means_flatten(:,ind);

	%90
	ti = [9 16  41 49  70 79  102 110  133 142  164 172  196 204  226 234]+0.5;
	%346
	ti = [9 17  41 49  71 79  103 111  134 141  164 173  196 204  226 234]+0.5;
	c0 = ones(1, 31*8);

	maxiter = 50000;
	lastcost = 1e10;
	currcost = lastcost;
	for iter=1:maxiter
		g = generate_alphas(tauguess, 8, ti);
		As = yactual.' / [c0; g];

		%ftau = @(tau) ftrial([tau, As]);
		ftau = @(tau) (As * [c0; generate_alphas(tau, 8, ti)]).'; % expansion of the above so we don't rely on ftrial which changes
		fftau = @(tau, xdata) (As * [c0; generate_alphas(tau, 8, ti)]).'; % expansion of the above so we don't rely on ftrial which changes
		%[tauguess, yfit, cost] = lsqsearch(ftau, tauguess, yactual);
		[tauguess, cost] = lsqcurvefit(fftau, tauguess, [], yactual);
		currcost = cost;
		if lastcost<currcost
			warning('whaaaaaat');
		elseif lastcost-currcost < 0.001
			break
		end
		lastcost = cost;
	end
	currcost
	g = generate_alphas(tauguess, 8, ti);
	As = yactual.' / [c0; g];
	yfit = (As * [c0; g]).';
	sum((yfit-yactual).^2)

	coeffs16exp_twostep(ind,:) = [tauguess, As];
	fit16exp_twostep(:,ind) = yfit;

end
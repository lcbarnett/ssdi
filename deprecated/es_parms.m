function esrule = es_parms(rule,dim)

% Set 1+1 evolutionary strategy parameters

if rule < eps  % special case (ignore gain): esrule = 1/5, gain = (5/4)*log(3/2) = 0.5068...
	ifac = 3/2;
	nfac = realpow(3/2,-1/4);
else
	gain = 1/sqrt(dim+1);
	ifac = exp((1-rule)*gain);
	nfac = exp(-rule*gain);
end
esrule = [ifac,nfac];

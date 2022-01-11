function gp_opthist(dhist,niters,titles,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

nhists = length(dhist);
for h = 1:nhists
	nruns(h) = length(dhist{h});
	ddmax{h} = max(cell2mat(dhist{h}));
	ddmin{h} = min(cell2mat(dhist{h}));
end
nnzhists = nnz(niters);

[~,gpname] = fileparts([gpstem '.xxx']); % hack to get fileparts to behave itself
for h = 1:nhists
	gp_write(sprintf('%s_%d',gpstem,h),dhist{h});
end
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
for h = 1:nhists
	fprintf(gp,'datfile_%d = "%s_%d.dat"\n\n',h,gpname,h);
end
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n\n');
fprintf(gp,'set multiplot title "%s" layout 2,%d\n\n',gptitle,nnzhists);

fprintf(gp,'set ylabel "DD" norot\n\n');
for h = 1:nhists
	if niters(h) == 0, continue; end
	fprintf(gp,'set title "%s - dynamical dependence"\n',titles{h});
	fprintf(gp,'set xr [1:%g]\n',niters(h));
	fprintf(gp,'set yr [0:%g]\n',1.05*ddmax{h}(2));
	fprintf(gp,'plot \\\n');
	for k = 1:nruns(h)
		fprintf(gp,'datfile_%d i %d using 1:2 w lines not, \\\n',h,k-1);
	end
	fprintf(gp,'NaN not\n\n');
end

fprintf(gp,'set ylabel "$\\\\sigma$" norot\n');
fprintf(gp,'set logs y\n\n');
for h = 1:nhists
	if niters(h) == 0, continue; end
	fprintf(gp,'set title "%s - step size"\n',titles{h});
	fprintf(gp,'set xr [1:%g]\n',niters(h));
	fprintf(gp,'set yr [*:%g]\n',1.05*ddmax{h}(3));
	fprintf(gp,'plot \\\n');
	for k = 1:nruns(h)
		fprintf(gp,'datfile_%d i %d using 1:3 w lines not, \\\n',h,k-1);
	end
	fprintf(gp,'NaN not\n\n');
end

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

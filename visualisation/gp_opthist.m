function gp_opthist(dhist,niters,havegrad,titles,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

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
fprintf(gp,'set lmargin 12\n');
fprintf(gp,'set rmargin 6\n');
fprintf(gp,'set grid\n\n');
fprintf(gp,'set multiplot title "%s" layout 3,%d\n\n',gptitle,nnzhists);

fprintf(gp,'set ylabel "DD" norot\n\n');
for h = 1:nhists
	if niters(h)> 0
		fprintf(gp,'set title "%s - dynamical dependence"\n',titles{h});
		fprintf(gp,'set xr [1:%g]\n',niters(h));
		fprintf(gp,'set yr [%g:%g]\n',0.95*ddmin{h}(1),1.05*ddmax{h}(1));
		fprintf(gp,'plot \\\n');
		for k = 1:nruns(h)
			fprintf(gp,'datfile_%d i %d using 1 w lines not, \\\n',h,k-1);
		end
		fprintf(gp,'NaN not\n\n');
	else
		fprintf(gp,'set multiplot next\n\n');
	end
end

fprintf(gp,'set ylabel "$\\\\sigma$" norot\n');
fprintf(gp,'set logs y\n\n');
for h = 1:nhists
	if niters(h) > 0
		fprintf(gp,'set title "%s - step size"\n',titles{h});
		fprintf(gp,'set xr [1:%g]\n',niters(h));
		fprintf(gp,'set yr [*:%g]\n',1.05*ddmax{h}(2));
		fprintf(gp,'plot \\\n');
		for k = 1:nruns(h)
			fprintf(gp,'datfile_%d i %d using 2 w lines not, \\\n',h,k-1);
		end
		fprintf(gp,'NaN not\n\n');
	else
		fprintf(gp,'set multiplot next\n\n');
	end
end

fprintf(gp,'set ylabel "$\\\\nabla$" norot\n');
fprintf(gp,'set logs y\n\n');
for h = 1:nhists
	if niters(h) > 0 && havegrad(h)
		fprintf(gp,'set title "%s - gradient"\n',titles{h});
		fprintf(gp,'set xr [1:%g]\n',niters(h));
		fprintf(gp,'set yr [*:%g]\n',1.05*ddmax{h}(2));
		fprintf(gp,'plot \\\n');
		for k = 1:nruns(h)
			fprintf(gp,'datfile_%d i %d using 3 w lines not, \\\n',h,k-1);
		end
		fprintf(gp,'NaN not\n\n');
	else
		fprintf(gp,'set multiplot next\n\n');
	end
end

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

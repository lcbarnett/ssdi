function gp_localopt(dopt,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot (local) optimum dynamical dependencies at all scales

[nruns,n1] = size(dopt);
domax = max(dopt(:));
precstr = cell(1,nruns);
precstr{1} = ' %4d';
for k = 2:nruns+1
	precstr{k} = ' %24.16f';
end

[~,gpname] = fileparts([gpstem '.xxx']); % hack to get fileparts to behave itself
gp_write(gpstem,[1:n1;dopt]',precstr);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set title "%s"\n',gptitle);
fprintf(gp,'unset key\n');
fprintf(gp,'set xr[-0.5:%g]\n',n1-0.5);
fprintf(gp,'set yr [0:%g]\n',1.05*domax);
fprintf(gp,'set style data histogram\n');
fprintf(gp,'set style histogram cluster gap 1\n');
fprintf(gp,'# set style fill solid 0.3 border lt -1\n');
fprintf(gp,'set style fill solid 0.4 noborder\n');
fprintf(gp,'set boxwidth 1\n');
fprintf(gp,'set xlabel "scale (m)"\n');
fprintf(gp,'set ylabel "optimal dynamical dependence"\n');
fprintf(gp,'set grid y\n');
fprintf(gp,'x = -0.5\n');
for k = 1:n1+1
	fprintf(gp,'set arrow from first x,graph 0 to first x,graph 1 nohead lt -1 lw 2; x = x+1\n');
end
fprintf(gp,'plot \\\n');
fprintf(gp,'datfile using 2:xtic(1) lc "#4D74C4" not, \\\n');
for k = 2:nruns
	fprintf(gp,'datfile using %d lc "#4D74C4" not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
gp_close(gp,gpstem,gpterm,gpplot);

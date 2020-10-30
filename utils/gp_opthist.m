function gp_opthist(dd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

[niters,nruns] = size(dd);
ii = (1:niters)';
ddmax = max(dd(:));

[~,gpname] = fileparts(gpstem);
gp_write(gpstem,[ii dd]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set title "%s"\n',gptitle);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xr [1:%g]\n',niters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddmax);
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set ylabel "dynamical dependence"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n');
fprintf(gp,'plot \\\n');
for k = 1:nruns
	fprintf(gp,'datfile using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
gp_close(gp,gpstem,gpterm,gpplot);

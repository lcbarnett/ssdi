function gp_opthist(ddp,npiters,ddo,noiters,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

npruns = length(ddp);
ddpmax = max(cell2mat(ddp));
ddpmin = min(cell2mat(ddp));

noruns = length(ddo);
ddomax = max(cell2mat(ddo));
ddomin = min(cell2mat(ddo));

[~,gpname] = fileparts(gpstem);
gp_write([gpstem '_p'],ddp);
gp_write([gpstem '_o'],ddo);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfilep = "%s.dat"\n\n',[gpname '_p']);
fprintf(gp,'datfileo = "%s.dat"\n\n',[gpname '_o']);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n\n');
fprintf(gp,'set multiplot title "%s" layout 2,2\n\n',gptitle);

fprintf(gp,'set ylabel "DD" norot\n\n');

fprintf(gp,'set title "Pre-optimisation - dynamical dependence"\n');
fprintf(gp,'set xr [1:%g]\n',npiters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddpmax(2));
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfilep i %d using 1:2 w lines not, \\\n',k-1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set title "Optimisation - dynamical dependence"\n');
fprintf(gp,'set xr [1:%g]\n',noiters);
fprintf(gp,'set yr [%g:%g]\n',0.95*ddomin(2),1.05*ddomax(2));
fprintf(gp,'plot \\\n');
for k = 1:noruns
	fprintf(gp,'datfileo i %d using 1:2 w lines not, \\\n',k-1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set ylabel "$\\\\sigma$" norot\n');
fprintf(gp,'set logs y\n\n');

fprintf(gp,'set title "Pre-optimisation - step size"\n');
fprintf(gp,'set xr [1:%g]\n',npiters);
fprintf(gp,'set yr [*:%g]\n',1.05*ddpmax(3));
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfilep i %d using 1:3 w lines not, \\\n',k-1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set title "Optimisation - step size"\n');
fprintf(gp,'set xr [1:%g]\n',noiters);
fprintf(gp,'set yr [*:%g]\n',1.05*ddomax(3));
fprintf(gp,'plot \\\n');
for k = 1:noruns
	fprintf(gp,'datfileo i %d using 1:3 w lines not, \\\n',k-1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

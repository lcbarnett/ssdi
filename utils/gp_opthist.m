function gp_opthist(ddp,ddn,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

[piters,npruns] = size(ddp);
ip = (1:piters)';
ddpmax = max(ddp(:));

[niters,nnruns] = size(ddn);
in = (1:niters)';
ddnmax = max(ddn(:));

[~,gpname] = fileparts(gpstem);
gp_write(gpstem,{[ip ddp],[in,ddn]});
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set ylabel "dynamical dependence"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n');
fprintf(gp,'set multiplot title "%s" layout 2,1\n',gptitle);
fprintf(gp,'set title "Pre-optimisation"\n');
fprintf(gp,'set xr [1:%g]\n',piters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddpmax);
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfile i 0 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
fprintf(gp,'set title "Optimisation"\n');
fprintf(gp,'set xr [1:%g]\n',niters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddnmax);
fprintf(gp,'plot \\\n');
for k = 1:nnruns
	fprintf(gp,'datfile i 1 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

function gp_opthist(ddp,dds,ddn,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

[piters,npruns] = size(ddp);
ip = (1:piters)';
ddpmax = max(ddp(:));
ddpmin = min(ddp(:));

[siters,nsruns] = size(dds);
is = (1:siters)';
ddsmax = max(dds(:));
ddsmin = min(dds(:));

[niters,nnruns] = size(ddn);
in = (1:niters)';
ddnmax = max(ddn(:));
ddnmin = min(ddn(:));

[~,gpname] = fileparts(gpstem);
gp_write(gpstem,{[ip ddp],[is dds],[in,ddn]});
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set ylabel "dynamical dependence"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n');
fprintf(gp,'set multiplot title "%s" layout 3,1\n',gptitle);
fprintf(gp,'set title "Pre-optimisation 1"\n');
fprintf(gp,'set xr [1:%g]\n',piters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddpmax);
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfile i 0 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
fprintf(gp,'set title "Pre-optimisation 2"\n');
fprintf(gp,'set xr [1:%g]\n',siters);
fprintf(gp,'set yr [%g:%g]\n',0.95*ddsmin,1.05*ddsmax);
fprintf(gp,'plot \\\n');
for k = 1:nsruns
	fprintf(gp,'datfile i 1 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
fprintf(gp,'set title "Optimisation"\n');
fprintf(gp,'set xr [1:%g]\n',niters);
fprintf(gp,'set yr [%g:%g]\n',0.95*ddnmin,1.05*ddnmax);
fprintf(gp,'plot \\\n');
for k = 1:nnruns
	fprintf(gp,'datfile i 2 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n');
fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

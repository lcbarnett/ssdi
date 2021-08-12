function gp_opthist(ddp,ddo,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot optimisation histories

[npiters,~,npruns] = size(ddp);
piters  = (1:npiters)';
ddp1    = squeeze(ddp(:,1,:));
ddp2    = squeeze(ddp(:,2,:));
ddp1max = max(ddp1(:));
ddp1min = min(ddp1(:));
ddp2max = max(ddp2(:));
ddp2min = min(ddp2(:));

[noiters,~,nsruns] = size(ddo);
oiters  = (1:noiters)';
ddo1    = squeeze(ddo(:,1,:));
ddo2    = squeeze(ddo(:,2,:));
ddo1max = max(ddo1(:));
ddo1min = min(ddo1(:));
ddo2max = max(ddo2(:));
ddo2min = min(ddo2(:));

[~,gpname] = fileparts(gpstem);
gp_write(gpstem,{[piters ddp1 ddp2],[oiters ddo1,ddo2]});
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n\n');
fprintf(gp,'set multiplot title "%s" layout 2,2\n\n',gptitle);

fprintf(gp,'set ylabel "DD" norot\n\n');

fprintf(gp,'set title "Pre-optimisation - dynamical dependence"\n');
fprintf(gp,'set xr [1:%g]\n',npiters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddp1max);
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfile i 0 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set title "Optimisation - dynamical dependence"\n');
fprintf(gp,'set xr [1:%g]\n',noiters);
fprintf(gp,'set yr [%g:%g]\n',0.95*ddo1min,1.05*ddo1max);
fprintf(gp,'plot \\\n');
for k = 1:nsruns
	fprintf(gp,'datfile i 1 using 1:%d w lines not, \\\n',k+1);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set ylabel "$\\\\sigma$" norot\n');
fprintf(gp,'set logs y\n\n');

fprintf(gp,'set title "Pre-optimisation - step size"\n');
fprintf(gp,'set xr [1:%g]\n',npiters);
fprintf(gp,'set yr [*:%g]\n',1.05*ddp2max);
fprintf(gp,'plot \\\n');
for k = 1:npruns
	fprintf(gp,'datfile i 0 using 1:%d w lines not, \\\n',k+1+npruns);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'set title "Optimisation - step size"\n');
fprintf(gp,'set xr [1:%g]\n',noiters);
fprintf(gp,'set yr [*:%g]\n',1.05*ddo2max);
fprintf(gp,'plot \\\n');
for k = 1:nsruns
	fprintf(gp,'datfile i 1 using 1:%d w lines not, \\\n',k+1+nsruns);
end
fprintf(gp,'NaN not\n\n');

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

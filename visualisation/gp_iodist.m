function gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot inter-optima subspace distances

nruns = size(Loptd,1);

[~,gpname] = fileparts([gpstem '.xxx']); % hack to get fileparts to behave itself
gp_write(gpstem,[0 1:nruns; (1:nruns)' Loptd]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set title "%s"\n',gptitle);
fprintf(gp,'set size square\n');
fprintf(gp,'unset key\n');
fprintf(gp,'unset xlabel\n');
fprintf(gp,'unset ylabel\n');
fprintf(gp,'set xr[0.5:%d+0.5]\n',nruns);
fprintf(gp,'set yr[0.5:%d+0.5]\n',nruns);
fprintf(gp,'set cbrange[0:1]\n');
fprintf(gp,'set palette defined ( 0 "#ffffff", 0.1 "#4D74C4", 1 "#000000" )\n');
fprintf(gp,'plot datfile matrix with image not\n');
gp_close(gp,gpstem,gpterm,gpplot);

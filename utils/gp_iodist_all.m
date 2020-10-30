function gp_iodist_all(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot inter-optima subspace distances

[nruns,~,n1] = size(Loptd);

precstr = ' %24.16f';
ncols = 2;
nrows = ceil(n1/ncols);

[~,gpname] = fileparts(gpstem);
gp = gp_open(gpstem,gpterm,[Inf,0.5],gpfsize);
fprintf(gp,'set size square\n');
fprintf(gp,'unset key\n');
fprintf(gp,'unset xlabel\n');
fprintf(gp,'unset ylabel\n');
fprintf(gp,'set xr[0.5:%d+0.5]\n',nruns);
fprintf(gp,'set yr[0.5:%d+0.5]\n',nruns);
fprintf(gp,'set cbrange[0:1]\n');
fprintf(gp,'set palette defined ( 0 "#ffffff", 0.6 "#4D74C4", 1 "#000000" )\n');
fprintf(gp,'set multiplot title "%s" layout %d,%d\n\n',gptitle,nrows,ncols);

for m = 1:n1
	fprintf(gp,'set title "scale = %d"\n',m);
	fprintf(gp,'plot "-" matrix with image not\n');
	fprintf(gp,precstr,0);
	for c = 1:nruns
		fprintf(gp,precstr,c);
	end
	fprintf(gp,'\n');
	for r = 1:nruns
		fprintf(gp,precstr,r);
		for c = 1:nruns
			fprintf(gp,precstr,Loptd(r,c,m));
		end
		fprintf(gp,'\n');
	end
	fprintf(gp,'EOD\n\n');
end

fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);

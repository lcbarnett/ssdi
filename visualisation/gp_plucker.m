function gp_plucker(Loptx,n,m,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot)

% Plot distances from subspaces spanned by all combinatios of m axes

c = nchoosek(1:n,m);
nc = size(c,1);

for k = 1:nc
	if Loptx(k) > 0.1
		fprintf('%3d : %s : %6.4f\n',k,num2str(c(k,:),' %d'),Loptx(k));
	end
end

[~,gpname] = fileparts(gpstem);
gp_write(gpstem,[(1:nc)' Loptx]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set title "%s"\n',gptitle);
fprintf(gp,'unset key\n');
fprintf(gp,'set xr[-0.5:%g]\n',nc-0.5);
fprintf(gp,'set yr [0:1.05]\n',1.05);
fprintf(gp,'set style data histogram\n');
fprintf(gp,'set style fill solid 1 noborder\n');
fprintf(gp,'set boxwidth 1\n');
fprintf(gp,'set xlabel "subspace"\n');
fprintf(gp,'set ylabel "normalised angle"\n');
fprintf(gp,'set grid\n');
fprintf(gp,'plot datfile u 1:2 w boxes lc "#4D74C4" not\n');
fprintf(gp,'\n');
for k = 1:nc
	if Loptx(k) > 0.1
		fprintf(gp,'# %3d : %s : %6.4f\n',k,num2str(c(k,:),' %d'),Loptx(k));
	end
end
fprintf(gp,'\n');
gp_close(gp,gpstem,gpterm,gpplot);

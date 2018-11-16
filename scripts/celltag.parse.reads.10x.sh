#!/usr/bin/gawk -f

BEGIN {
	FS="\t";
	OFS="\t";
	sep=":"
	xcregex = "CB:Z:[ACTGN0-9]+"
	xmregex = "UB:Z:[ACTGN0-9]+"
	rgregex = "GN:Z:[A-Za-z0-9]+"
	#tagregex = "CCGGT([ACTG]{8})GAATTC";
	print "Read.ID", "Read.Seq", "Cell.BC", "UMI", "Cell.Tag", "Gene"
}
{
	match($0, xcregex);
	xc = substr($0, RSTART, RLENGTH);
	xc = split(xc, CBC, sep);
	xc = CBC[3];
	if(!xc){
		match($0, "CR:Z:[ACTGN0-9]+");
		xc = substr($0, RSTART, RLENGTH);
		xc = split(xc, CBC, sep)
		xc = CBC[3]
	}

	match($0, xmregex);
	xm = substr($0, RSTART, RLENGTH);
	xm = split(xm, UMI, sep);
	xm = UMI[3]
	match($0, rgregex);
	rg = substr($0, RSTART, RLENGTH);
	rg = split(rg, GENE, sep);
	rg = GENE[3]
	if(rg == "" || rg == "\t"){
		rg = "BAB"
	}
	match($10, tagregex, celltag);
	#if( xc && xm && celltag[1] ){
		print $1, $10, xc, xm, celltag[1], rg;
	#}
}

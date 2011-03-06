#!/bin/bash

echo "Please remove all the lines above #!/usr/bin/awk -f"
echo "This is just because CVS is configured to not accept a script beginning with #!/usr/bin/awk -f"
echo "and #!/usr/bin/env awk -f doesn't do the job"

exit 0

#!/usr/bin/awk -f

BEGIN {
	# Check arguments
	if (ARGC < 1) {
		print "Argument error! Specify input file.";
		exit(1);
	}

	# Remove argument from ARGV and add implied "-" (stdin)
	if (ARGC == 0)
		ARGV[0] = "-";
	else
		ARGV[0] = "";

    i = 0;
}

/^#/ { next; }

{
    ptbin = sprintf("%4.1f -- %4.1f", $1,$2);
    etabin = sprintf("%4.2f -- %4.2f",$3,$4);
    ptbins[ptbin]=i;
    etabins[etabin]=i;
    eff[ptbin,etabin] = $7;
    efferr[ptbin,etabin] = $8;
    i++
}

END {
    # make an reversed map and sort it
    for (a in ptbins) { ptbinsRev[ptbins[a]]=a; }
    nPtbins = asort(ptbinsRev);

    for (a in etabins) { etabinsRev[etabins[a]]=a; }
    nEtabins = asort(etabinsRev);

    printf "\\begin{tabular}{c|";
    for (i=0; i!=length(etabins); i++)
	printf "c";
    printf "}\n";
    print "\\dbline";
    printf "\\pt & \\multicolumn{%d}{c}{$\\eta$} \\\\ \n", length(etabins);
    printf "GeV/c";
    for (i = 1; i<= nEtabins; i++) {
	printf " & %s", etabinsRev[i];
    }
    printf "\\\\\n";
    print "\\sgline";
    for (i = 1; i<= nPtbins; i++) {
	ptbin = ptbinsRev[i];
	printf "%s ", ptbin;
	for (j = 1; j<=nEtabins; j++) {
	    etabin = etabinsRev[j];
	    val = eff[ptbin,etabin];
	    err = efferr[ptbin,etabin];
	    if (val < 0)
		printf " &  ---  ";
	    else
		printf " & \\ensuremath{%5.3f \\pm %5.3f}", val, err;
	}
	printf "\\\\\n";
    }
    print "\\dbline";
    print "\\end{tabular}";
}


mmseqs createdb query.fasta queryDB
mmseqs createdb target.fasta targetDB
mmseqs search queryDB targetDB resultDB tmp -a -s 11 -e inf --min-ungapped-score 0 --zdrop 100000 --diag-score 0 --k-score 25
mmseqs convertalis queryDB targetDB resultDB resultDB.m8

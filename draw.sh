supnet -n "$1" -d | dot -Tpdf > p.pdf && evince p.pdf

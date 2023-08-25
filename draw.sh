supnet_dtcache -n "$1" --dot | dot -Tpdf > p.pdf && evince p.pdf

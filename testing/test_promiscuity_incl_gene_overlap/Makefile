test:
	../../util/promiscuity_filter.pl  --fusion_preds star-fusion.fusion_candidates.preliminary --genome_lib_dir ${CTAT_GENOME_LIB} --max_promiscuity 2
	diff star-fusion.fusion_candidates.preliminary.post_promisc_filter __prom2_out

	../../util/promiscuity_filter.pl  --fusion_preds star-fusion.fusion_candidates.preliminary --genome_lib_dir ${CTAT_GENOME_LIB} --max_promiscuity 3
	diff star-fusion.fusion_candidates.preliminary.post_promisc_filter __prom3_out

clean:
	./cleanMe.sh


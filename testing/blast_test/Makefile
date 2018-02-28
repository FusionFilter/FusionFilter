test:
	./runMe.sh
	sort fusion_preds.txt.post_blast_filter > fusion_preds.txt.post_blast_filter.sorted
	diff fusion_preds.txt.post_blast_filter.sorted fusion_preds.txt.expected_result


clean:
	./cleanMe.sh


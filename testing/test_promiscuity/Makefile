

test:
	./runMe.sh
	sort fusion_preds.txt.post_promisc_filter > fusion_preds.txt.post_promisc_filter.sorted
	diff fusion_preds.txt.post_promisc_filter.sorted fusion_preds.txt.expected_result

clean:
	./cleanMe.sh

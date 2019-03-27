DIRS = blast_test  \
       test_promiscuity \
       test_promiscuity_incl_gene_overlap \
       test_paralogs

test:
	@for i in $(DIRS); do \
	echo "Running test in $$i..."; \
	(cd $$i; $(MAKE) test) || exit $$?; done

clean:
	@for i in $(DIRS); do \
	echo "Running test in $$i..."; \
	(cd $$i; $(MAKE) clean) || exit $$?; done


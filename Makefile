.PHONY: test style run rrv

test:
	python3 poisson_cat.py

style:
	black run.py -l 79

# Shows an example of how to use this in conjunction with rankratioviz
run:
	./run.py --table /projects/emp/redbiom/emp_deblur_150bp.release1.biom \
			 --metadata /projects/emp/redbiom/emp_deblur_150bp.release1.map.tsv \
			 --category empo_1 \
			 --output-path diffs.tsv \
			 --filter-control

rrv:
	rm -rf emp_poisson_rrv_out/
	rankratioviz --table /projects/emp/redbiom/emp_deblur_150bp.release1.biom \
			 	 --ranks diffs.tsv \
				 --sample-metadata /projects/emp/redbiom/emp_deblur_150bp.release1.map.tsv \
				 --output-dir emp_poisson_rrv_out/

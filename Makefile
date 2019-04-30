.PHONY: test style run_emp rrv_emp

test:
	python3 poisson_cat.py

# Note: this requires that "black" is installed.
style:
	black run.py -l 79

# Shows an example of how to use this in conjunction with rankratioviz
run_emp:
	./run.py \
		--table ../EMP/emp_deblur_150bp.release1.biom \
		--metadata ../EMP/emp_deblur_150bp.release1.map.tsv \
		--category "empo_1" \
		--output-path diffs.tsv \
		--filter-category-value "Control"

rrv_emp:
	rm -rf emp_poisson_rrv_out/
	rankratioviz \
		--table /projects/emp/redbiom/emp_deblur_150bp.release1.biom \
		--ranks diffs.tsv \
		--sample-metadata /projects/emp/redbiom/emp_deblur_150bp.release1.map.tsv \
		--output-dir emp_poisson_rrv_out/

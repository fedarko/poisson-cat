# poisson-cat

Closed-form poisson differential computation.

All of the code in `poisson_cat.py` was written by Jamie Morton (@mortonjt)
and Cameron Martino (@cameronmartino), who also worked out the accompanying
math.

## Running this

This code is dependent on the following python packages:
    - click
    - biom-format
    - pandas
    - numpy
    - scipy
    - scikit-learn

```
Usage: run.py [OPTIONS]

Options:
  -t, --table TEXT                BIOM table with count data  [required]
  -m, --metadata TEXT             Sample metadata file  [required]
  -c, --category TEXT             Metadata category of interest  [required]
  -r, --reference-category TEXT   Reference metadata category of interest; if
                                  not specified, the first category will be
                                  picked
  -o, --output-path TEXT          Output filepath to which differentials TSV
                                  will be written  [required]
  -f, --filter-control / --no-filter-control
                                  If passed, this will filter out all samples
                                  with a category value of "Control"; will
                                  also then filter out all "empty" features
  --help                          Show this message and exit.
```

The output differentials should be displayable in
[rankratioviz](https://github.com/fedarko/rankratioviz).

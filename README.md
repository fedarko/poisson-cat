# poisson-cat
[![Build Status](https://travis-ci.org/fedarko/poisson-cat.svg?branch=master)](https://travis-ci.org/fedarko/poisson-cat)

Closed-form poisson differential computation.

**All of the code in `poisson_cat.py` was written by Jamie Morton
([@mortonjt](https://github.com/mortonjt))
and Cameron Martino ([@cameronmartino](https://github.com/cameronmartino)),
who also worked out the accompanying math.** See
[here](https://gist.github.com/mortonjt/0bb8d0565fc6c02b48e91524e816f112) for
the original gist from Jamie, and
[here](https://gist.github.com/cameronmartino/8d89b73c2dcd749992127ad5a8d284e2)
for Cameron's version of it.

## Installation

```bash
pip install click biom-format pandas numpy scipy scikit-learn
```

## Running this

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

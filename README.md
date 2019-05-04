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
  -c, --category TEXT             Metadata category of interest; currently
                                  only binary categories (i.e. those
                                  containing only two unique values) are
                                  supported  [required]
  -r, --reference-category TEXT   Reference metadata category of interest; if
                                  not specified, the first category will be
                                  picked
  -o, --output-path TEXT          Output filepath to which differentials TSV
                                  will be written  [required]
  -f, --filter-category-value TEXT
                                  If passed, this will filter out all samples
                                  with a -c category value of this string.
                                  This will also afterwards filter out all
                                  "empty" features. This is useful if you have
                                  a category with three possible values that
                                  you'd like to make into a binary category,
                                  so that it can be used here.
  --help                          Show this message and exit.
```

## Visualizing the output

The output differentials can be displayed in
[rankratioviz](https://github.com/fedarko/rankratioviz). If you're using
rankratioviz through QIIME 2, you can import the output differentials here as a
`FeatureData[Differential]` artifact, and visualize them using
`qiime rankratioviz supervised-rank-plot`; if you're using rankratioviz
outside of QIIME 2, you should just be able to pass the differentials directly
to rankratioviz.

(Note that you might want to use the `-x / --extreme-feature-count` option in
rankratioviz if your BIOM table has a lot of entries.)

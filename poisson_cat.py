import numpy as np
import pandas as pd
from biom import Table
from scipy.stats import pearsonr
from sklearn.utils import check_random_state

# This is the main function


def poisson_cat(table, metadata, category, ref=None):
    """
    This method calculates Poisson differential abundance
    for binary categorical covariates.

    Parameters
    ----------
    table : biom.Table
        Sparse matrix of counts.
    metadata : pd.DataFrame
        Sample metadata file.
    category : str
        Metadata category of interest.
    ref : str
        Reference category.  If not specified, the first
        category will be picked.

    Returns
    -------
    differential : pd.DataFrame
        Differentials of all of the variables of interest.

    Raises
    ------
    ValueError
        The sample names are inconsistent
        between `metadata` and `table`.
        Please ensure they are matched.
    ValueError
        Less than two catagories given.
        Please only provide binary catagories.
    ValueError
        More than two catagories given.
        Please only provide binary catagories.

    TODO
    ----
    Extend this to multiple metadata categories, so that multiple columns
    can be analyzed simultaneously
    """

    # test that metadata and table match
    if len(set(table.ids(axis='sample')) & set(
            metadata.index)) != len(metadata.index):
        # not matched before using
        # return error, needs to be
        # matched externally
        raise ValueError('The sample names are inconsistent',
                         'between `metadata` and `table`.',
                         'Please ensure they are matched.')
    # get the set of sub-catagories
    cats = metadata[category].unique()
    # test that category is binary
    if cats.shape[0] < 2:
        raise ValueError('Less than two catagories given.',
                         'Please only provide binary catagories.')
    if cats.shape[0] > 2:
        raise ValueError('More than two catagories given.',
                         'Please only provide binary catagories.')
    # if cats not [-1,1] convert it for calculation
    cat_map = {y: x for x, y in zip([-1, 1], cats)}
    metadata[category] = [cat_map[x] for x in metadata[category]]
    # if no reference provided use first
    if ref is None:
        idx = metadata[category] == cat_map[cats[0]]
    # if provided use that one
    else:
        idx = metadata[category] == cat_map[ref]
    # calculate fast-poisson
    x = idx.astype(np.int64)  # sample covariates
    y = table.matrix_data  # table counts
    A = y @ ~x  # reference microbe counts
    # total counts by microbe
    B = table.sum(axis='observation')
    # total counts by sample
    N = table.sum(axis='sample')
    # total counts in reference samples
    C = N[idx].sum()
    # total counts in non-ref samples
    D = N[~idx].sum()
    # calculate the differential - ensure non-negative (breaks log)
    diff = pd.Series(np.log(-1 * ((A * B * C) / (B * B * D - A))),
                     table.ids(axis='observation'))
    return diff


# This is for testing
def random_block_table(reps, n_species,
                       species_mean=0,
                       species_var=1.,
                       effect_size=1,
                       library_size=10000,
                       microbe_total=100000, microbe_kappa=0.3,
                       microbe_tau=0.1, sigma=0.5, seed=None):
    """ Differential abundance analysis benchmarks.

    The simulation here consists of 3 parts

    Step 1: generate class probabilities using logistic distribution
    Step 2: generate coefficients from normal distributions
    Step 3: generate counts from species distributions

    Parameters
    ----------
    reps : int
        Number of replicate samples per test.
    n_species : int
        Number of species.
    species_loc : float
        Mean of the species prior.
    species_variance : float
        Variance of species log-fold differences
    effect_size : int
        The effect size difference between the feature abundances.
    n_contaminants : int
       Number of contaminant species.
    sigma: float
        Logistic error variance for class probabilities
    library_size : np.array
        A vector specifying the library sizes per sample.
    template : np.array
        A vector specifying feature abundances or relative proportions.

    Returns
    -------
    generator of
        pd.DataFrame
           Ground truth tables.
        pd.DataFrame
           Metadata group categories, n_diff and effect_size
        pd.Series
           Species actually differentially abundant.
    """
    state = check_random_state(seed)

    n = reps * 2
    k = 2
    labels = np.array([-effect_size] * (n // 2) + [effect_size] * (n // 2))
    eps = np.random.logistic(loc=0, scale=sigma, size=n)
    class_probs = labels + eps

    X = np.hstack((np.ones((n, 1)), class_probs.reshape(-1, 1)))
    B = np.random.normal(
        loc=species_mean,
        scale=species_var,
        size=(
            k,
            n_species))

    # Helper functions
    # Convert microbial abundances to counts
    def to_counts_f(x):
        n = state.lognormal(np.log(library_size), microbe_tau)
        p = x / x.sum()
        return state.poisson(state.lognormal(np.log(n * p), microbe_kappa))

    o_ids = ['F%d' % i for i in range(n_species)]
    s_ids = ['S%d' % i for i in range(n)]

    abs_table = pd.DataFrame(np.exp(X @ B) * microbe_total,
                             index=s_ids,
                             columns=o_ids)

    rel_data = np.vstack(abs_table.apply(to_counts_f, axis=1))

    rel_table = pd.DataFrame(rel_data,
                             index=s_ids,
                             columns=o_ids)

    metadata = pd.DataFrame({'labels': labels})
    metadata['effect_size'] = effect_size
    metadata['microbe_total'] = microbe_total
    metadata['class_logits'] = class_probs
    metadata['intercept'] = 1
    metadata.index = s_ids

    ground_truth = pd.DataFrame({
        'intercept': B[0, :],
        'categorical': B[1, :]
    }, index=o_ids)

    return abs_table, rel_table, metadata, ground_truth


if __name__ == "__main__":

    import unittest

    class TestPoissonCat(unittest.TestCase):

        def setUp(self):
            reps = 50
            n_species = 200
            np.random.seed(2)
            self.res = random_block_table(
                reps, n_species,
                species_mean=0,
                species_var=1,
                microbe_kappa=0.7,
                microbe_tau=0.7,
                library_size=10000,
                microbe_total=100000,
                effect_size=1)
            abs_table, rel_table, metadata, ground_truth = self.res
            self.table = Table(rel_table.values.T,
                               rel_table.columns,
                               rel_table.index)

        def test_poisson_cat(self):
            abs_table, rel_table, metadata, ground_truth = self.res
            exp_diff = ground_truth['categorical']
            # test with no reference
            res_diff = poisson_cat(self.table, metadata,
                                   category='labels')
            r, p = pearsonr(res_diff, exp_diff)
            self.assertGreater(r**2, 0.5)
            self.assertLess(p, 0.01)
            # test with a reference
            res_diff_ref = poisson_cat(self.table, metadata,
                                       category='labels', ref=-1)
            r, p = pearsonr(res_diff_ref, exp_diff)
            self.assertGreater(r**2, 0.5)
            self.assertLess(p, 0.01)
            # ensure results same
            np.testing.assert_array_almost_equal(res_diff,
                                                 res_diff_ref)

        def test_poisson_cat_str(self):
            abs_table, rel_table, metadata, ground_truth = self.res
            exp_diff = ground_truth['categorical']
            mapflip = {x: y for x, y in zip([-1, 1], ['Sick', 'Healthy'])}
            metadata['labels'] = [mapflip[i] for i in metadata['labels']]
            # test with no reference
            res_diff = poisson_cat(self.table, metadata,
                                   category='labels')
            r, p = pearsonr(res_diff, exp_diff)
            self.assertGreater(r**2, 0.5)
            self.assertLess(p, 0.01)

        def test_poisson_cat_num(self):
            abs_table, rel_table, metadata, ground_truth = self.res
            exp_diff = ground_truth['categorical']
            metadata['labels'] += np.random.randint(20)
            # test with no reference
            res_diff = poisson_cat(self.table, metadata,
                                   category='labels')
            r, p = pearsonr(res_diff, exp_diff)
            self.assertGreater(r**2, 0.5)
            self.assertLess(p, 0.01)

        def test_poisson_multiclass(self):
            metadata = self.res[2]
            # test three catagories raises ValueError
            with self.assertRaises(ValueError):
                ridx = np.random.choice(metadata.index)
                metadata.loc[ridx, 'labels'] = 0.0
                poisson_cat(self.table, metadata,
                            category='labels')

        def test_poisson_oneclass(self):
            metadata = self.res[2]
            # test one category raises ValueError
            with self.assertRaises(ValueError):
                metadata = metadata.copy()
                metadata.loc[metadata.labels == -1, 'labels'] = 1.0
                poisson_cat(self.table, metadata,
                            category='labels')

        def test_poisson_notmatched(self):
            metadata = self.res[2]
            # test non-matched raises ValueError
            with self.assertRaises(ValueError):
                ridx = np.random.choice(metadata.index)
                metadata.drop(ridx, inplace=True)
                poisson_cat(self.table, metadata,
                            category='labels')

    unittest.main()

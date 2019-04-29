#!/usr/bin/env python3
import click
from biom import load_table
import pandas as pd
from poisson_cat import poisson_cat


# def get_row_indices_where_value_in_col(df, col_name, value):
#    row_indices = []
#    i = 0
#    for row_val in df[col_name]:
#        if row_val == value:
#            row_indices.append(i)
#        i += 1
#    return row_indices


@click.command()
@click.option(
    "-t", "--table", required=True, help="BIOM table with count data"
)
@click.option("-m", "--metadata", required=True, help="Sample metadata file")
@click.option(
    "-c", "--category", required=True, help="Metadata category of interest"
)
@click.option(
    "-r",
    "--reference-category",
    default=None,
    help=(
        "Reference metadata category of interest; if not "
        "specified, the first category will be picked"
    ),
)
@click.option(
    "-o",
    "--output-path",
    required=True,
    help="Output filepath to which differentials TSV will be written",
)
@click.option(
    "-f",
    "--filter-control/--no-filter-control",
    default=False,
    help=(
        "If passed, this will filter out all samples with a category value "
        'of "Control"'
    ),
)
def run_poisson_cat(
    table: str,
    metadata: str,
    category: str,
    reference_category: str,
    output_path: str,
    filter_control: bool,
) -> None:

    # table_df = load_table(table).to_dataframe()
    loaded_table = load_table(table)
    metadata_df = pd.read_csv(metadata, index_col=0, sep="\t")
    unique_cats = metadata_df[category].unique()
    if filter_control and unique_cats.shape[0] > 2:
        if "Control" in unique_cats:
            # Based on https://stackoverflow.com/a/18173074/1073
            print(
                "Number of samples pre-filtering: {}".format(
                    metadata_df.shape[0]
                )
            )
            control_row_idxs = metadata_df[
                metadata_df[category] == "Control"
            ].index
            metadata_df = metadata_df[metadata_df[category] != "Control"]
            loaded_table.filter(metadata_df.index)
            print(
                "Number of samples post-filtering: {}".format(
                    metadata_df.shape[0]
                )
            )
            print(
                "Number of features pre-filtering those with 0 counts: {}".format(
                    loaded_table.shape[0]
                )
            )
            # remove features in table with 0 counts (to eliminate any features
            # that were only present in now-filtered-out samples)
            loaded_table.remove_empty(axis="observation")
            print(
                "Number of features post-filtering those with 0 counts: {}".format(
                    loaded_table.shape[0]
                )
            )
            # TODO remove samples in table without a certain amount of reads
            # supporting them?

    print("Running poisson_cat...")
    diff = poisson_cat(loaded_table, metadata_df, category, reference_category)
    print("Done.")
    diff.to_csv(output_path, sep="\t", header=["Differential"], index_label="FeatureID")


if __name__ == "__main__":
    run_poisson_cat()

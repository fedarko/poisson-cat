#!/usr/bin/env python3
import click
from biom import load_table
import pandas as pd
from poisson_cat import poisson_cat


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
def run_poisson_cat(
    table: str,
    metadata: str,
    category: str,
    reference_category: str,
    output_path: str,
) -> None:
    loaded_table = load_table(table)
    metadata_df = pd.read_csv(metadata, index_col=0, sep="\t")
    print("Running poisson_cat...")
    diff = poisson_cat(loaded_table, metadata_df, category, reference_category)
    print("Done.")
    diff.to_csv(output_path, sep="\t")


if __name__ == "__main__":
    run_poisson_cat()

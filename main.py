from lib.deletion import DeletionSet, plot_deletion_set, read_deletion_sets_from_file
from lib.merge_priority_intervals import merge_deletion_set
from lib.cli import CliArgs
from lib.vcf import write_vcf_file
from typing import Callable
from concurrent.futures import ProcessPoolExecutor
import logging
import argparse
import time

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

type DeletionSetMergerFunction = Callable[[DeletionSet], DeletionSet]


def process_deletion_set(
    args: tuple[DeletionSet, DeletionSetMergerFunction, CliArgs],
) -> DeletionSet:
    deletion_set, deletion_set_merger, cli_args = args
    logger.info(
        f"Merging deletion set: {deletion_set.line_sub} ({deletion_set.chromosome})"
    )
    if cli_args.plot_diagrams:
        plot_deletion_set("unmerged_deletion_plots", deletion_set)
    merged = deletion_set_merger(deletion_set)
    if cli_args.plot_diagrams:
        plot_deletion_set("merged_deletion_plots", merged)
    return merged


def merge_deletion_sets(
    deletion_set_merger: DeletionSetMergerFunction, cli_args: CliArgs
) -> list[DeletionSet]:
    deletion_sets = list(read_deletion_sets_from_file())
    with ProcessPoolExecutor() as executor:
        merged_deletion_sets = list(
            executor.map(
                process_deletion_set,
                [(ds, deletion_set_merger, cli_args) for ds in deletion_sets],
            )
        )
    return merged_deletion_sets


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="zygofuser",
        description="Merges homozygous and hemizygous deletion lines, giving priority to homozygous calls.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        default=False,
        help="Plot deletion lines before and after merging",
    )
    args = parser.parse_args()
    cli_args = CliArgs(args.plot)

    tic = time.perf_counter()
    merged_deletion_sets = merge_deletion_sets(merge_deletion_set, cli_args)
    write_vcf_file(merged_deletion_sets)
    toc = time.perf_counter()
    logging.info(f"Done in {toc - tic} seconds")


if __name__ == "__main__":
    main()

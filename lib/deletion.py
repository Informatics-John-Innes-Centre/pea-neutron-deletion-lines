import os
import csv
from dataclasses import dataclass
import matplotlib.pyplot as plt
import matplotlib.patches as patches


@dataclass
class Deletion:
    start: int
    end: int
    homozygous: bool


class InvertedDeletionError(ValueError):
    def __init__(self, deletion: Deletion, chromosome: str, line_sub: str):
        self.deletion = deletion
        self.chromosome = chromosome
        self.line_sub = line_sub

        super().__init__(
            f"Deletion set for '{line_sub}' in chromosome '{chromosome}' contains inverted {'homozygous' if deletion.homozygous else 'hemizygous'} deletion of range {deletion.start} - {deletion.end}"
        )


@dataclass
class DeletionSet:
    deletions: list[Deletion]
    line_sub: str
    chromosome: str


def plot_deletion_set(folder: str, deletion_set: DeletionSet) -> None:
    os.makedirs(folder, exist_ok=True)
    _len = len(deletion_set.deletions)
    if _len == 0:
        return

    fig, ax = plt.subplots(figsize=(4, max(1.2, _len * 0.3 + 0.6)))
    for i, d in enumerate(deletion_set.deletions):
        ax.barh(
            i,
            d.end - d.start,
            left=d.start,
            height=0.6,
            color="#d62728" if not d.homozygous else "#1f77b4",
        )
    ax.set_yticks([])
    ax.set_ylim(-0.5, _len - 0.5)
    for i in range(_len - 1):
        ax.axhline(i + 0.5, color="#cccccc", linewidth=0.8, zorder=0)
    ax.set_xlabel("Position", fontsize=7)
    ax.set_title(
        f"{deletion_set.line_sub} - {deletion_set.chromosome}",
        fontsize=8,
        fontweight="bold",
    )
    ax.tick_params(axis="x", labelsize=6)
    ax.spines[["top", "right", "left"]].set_visible(False)
    hemi_patch = patches.Patch(color="#1f77b4", alpha=0.7, label="Homozygous")
    homo_patch = patches.Patch(color="#d62728", alpha=0.7, label="Hemizygous")
    ax.legend(
        handles=[hemi_patch, homo_patch],
        fontsize=6,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(
        f"{folder}/{deletion_set.line_sub}-{deletion_set.chromosome}.png", dpi=120
    )
    plt.close(fig)


def read_all_line_subs() -> list[str]:
    line_subs = []
    with open("data/metadata.csv") as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            line_subs.append(row[1])
    return sorted(set(line_subs))


def read_deletion_sets_from_file() -> list[DeletionSet]:
    # gather all our deletions
    deletion_sets: dict[tuple[str, str], DeletionSet] = {}
    with open("data/supplementary_table_5_excerpt.csv") as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # skip the header row
        for row in reader:
            # empty
            if len(row[5]) == 0:
                continue
            # get data for hemizygous
            if row[2]:
                line_sub = row[0].replace("/", "_")
                chromosome = row[2]
                start = int(row[3])
                end = int(row[4])
                homozygous = False

                if (line_sub, chromosome) in deletion_sets:
                    deletion_sets[(line_sub, chromosome)].deletions.append(
                        Deletion(start, end, homozygous)
                    )
                else:
                    deletion_sets[(line_sub, chromosome)] = DeletionSet(
                        [Deletion(start, end, homozygous)], line_sub, chromosome
                    )

            # get data for homozygous
            if row[8]:
                chromosome = row[8]
                start = int(row[9])
                end = int(row[10])
                line_sub = row[6].replace("/", "_")
                homozygous = True

                if (line_sub, chromosome) in deletion_sets:
                    deletion_sets[(line_sub, chromosome)].deletions.append(
                        Deletion(start, end, homozygous)
                    )
                else:
                    deletion_sets[(line_sub, chromosome)] = DeletionSet(
                        [Deletion(start, end, homozygous)], line_sub, chromosome
                    )
    return list(deletion_sets.values())

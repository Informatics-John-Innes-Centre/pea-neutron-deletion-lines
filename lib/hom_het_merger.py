"""Old algorithm for merging hemi/homozygous deletions"""

from lib.deletion import Deletion, DeletionSet, plot_deletion_set
from copy import deepcopy


def merge_overlapping_deletions_of_same_zygosity(
    deletions: list[Deletion],
) -> list[Deletion]:
    """
    Merge all overlapping deletions of the SAME zygosity - returns deletions merged in no particular order.
    """
    homoyzgous = [deletion for deletion in deletions if deletion.homozygous]
    hemizygous = [deletion for deletion in deletions if not deletion.homozygous]

    def merge(deletions: list[Deletion]) -> list[Deletion]:
        if not deletions:
            return []
        # Sort by start
        deletions.sort(key=lambda x: x.start)
        merged = [deletions[0]]
        for current in deletions[1:]:
            last = merged[-1]
            if current.start <= last.end:  # Overlap: extend if needed
                merged[-1].end = max(last.end, current.end)
            else:  # No overlap: add new interval
                merged.append(Deletion(current.start, current.end, current.homozygous))
        return merged

    return merge(homoyzgous) + merge(hemizygous)


def hom_het_merger(deletion_set: DeletionSet) -> DeletionSet:
    """
    Merge a set of deletions, prioritizing homozygous wherever it is present.
    First any overlapping homozygous to homozygous and hemizygous to hemizygous deletions will be merged,
    then the homozygous to hemizygous and vice versa will be merged.
    """
    # just take a copy for petes sake
    deletion_set = deepcopy(deletion_set)
    # merge any overlapping deletions of same zygosity
    deletion_set.deletions = merge_overlapping_deletions_of_same_zygosity(
        deletion_set.deletions
    )
    # sort by start pos, tie braker is zygosity - homozygous comes first
    deletion_set.deletions = sorted(
        deletion_set.deletions,
        key=lambda deletion: (deletion.start, not deletion.homozygous),
    )
    # merge deletions of opposite zygosity
    merged_deletions = [deletion_set.deletions[0]]
    for current_deletion in deletion_set.deletions[1:]:
        last_deletion = merged_deletions[-1]
        should_merge = current_deletion.homozygous != last_deletion.homozygous and (
            last_deletion.start <= current_deletion.end <= last_deletion.end
            or last_deletion.start <= current_deletion.start <= last_deletion.end
        )
        if should_merge:
            if not last_deletion.homozygous and current_deletion.homozygous:
                # Hemizygous First!
                # Hemizygous engulphs homozygous - we split it
                if (
                    last_deletion.start <= current_deletion.start
                    and last_deletion.end >= current_deletion.end
                ):
                    save_end = last_deletion.end
                    last_deletion.end = current_deletion.start
                    merged_deletions.append(current_deletion)
                    if current_deletion.end != save_end:
                        new_deletion = Deletion(current_deletion.end, save_end, False)
                        merged_deletions.append(new_deletion)
                    continue
                # hemizygous left overlap
                if (
                    last_deletion.start < current_deletion.start
                    and last_deletion.end > current_deletion.start
                    and last_deletion.end < current_deletion.end
                ):
                    last_deletion.end = current_deletion.start
                    merged_deletions.append(current_deletion)
                    continue
            elif last_deletion.homozygous and not current_deletion.homozygous:
                # Homozygous First!
                if (
                    last_deletion.start <= current_deletion.start
                    and last_deletion.end >= current_deletion.end
                ):
                    continue
                # homozygous left overlap hemizygous
                if (
                    last_deletion.start <= current_deletion.start
                    and last_deletion.end > current_deletion.start
                    and last_deletion.end < current_deletion.end
                ):
                    current_deletion.start = last_deletion.end
                    merged_deletions.append(current_deletion)
                    continue
            # Getting here is bad
            plot_deletion_set(
                "error",
                DeletionSet(
                    [last_deletion, current_deletion],
                    deletion_set.line_sub,
                    deletion_set.chromosome,
                ),
            )
            raise ValueError(
                f"Unhandled merging case {last_deletion}, {current_deletion}"
            )
        else:
            merged_deletions.append(current_deletion)
    return DeletionSet(merged_deletions, deletion_set.line_sub, deletion_set.chromosome)

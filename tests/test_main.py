import pytest
from lib.deletion import Deletion, DeletionSet, plot_deletion_set, InvertedDeletionError
from lib.merge_priority_intervals import merge_deletion_set


def both_forward_and_reverse_merge_helper(
    input: list[Deletion], expected_output: list[Deletion], name: str
) -> None:
    forward_deletion_set = DeletionSet(input, name + "-forward", "chr1")
    reversed_deletion_set = DeletionSet(input[::-1], name + "-reversed", "chr1")
    output_deletion_set = DeletionSet(expected_output, name, "chr1")
    plot_deletion_set("test_unmerged", forward_deletion_set)
    plot_deletion_set("test_unmerged", reversed_deletion_set)
    forward_merged = merge_deletion_set(forward_deletion_set)
    plot_deletion_set("test_merged", forward_merged)
    assert forward_merged.deletions == output_deletion_set.deletions
    reverse_merged = merge_deletion_set(reversed_deletion_set)
    plot_deletion_set("test_merged", reverse_merged)
    assert reverse_merged.deletions == output_deletion_set.deletions


def helper(input: list[Deletion], expected_output: list[Deletion], name: str) -> None:
    deletion_set = DeletionSet(input, name, "chr1")
    output_deletion_set = DeletionSet(expected_output, name, "chr1")
    plot_deletion_set("test_unmerged", deletion_set)
    merged = merge_deletion_set(deletion_set)
    plot_deletion_set("test_merged", merged)
    assert merged == output_deletion_set


def test_merge_homo_engulph_hemi() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 7, True), Deletion(2, 4, False)],
        [Deletion(0, 7, True)],
        "homo_engulph_hemi",
    )


def test_merge_homo_left_overlap() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(2, 7, False)],
        [Deletion(0, 5, True), Deletion(5, 7, False)],
        "homo_left_overlap_hemi",
    )


def test_merge_hemi_engulph_homo() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 7, False), Deletion(2, 4, True)],
        [Deletion(0, 2, False), Deletion(2, 4, True), Deletion(4, 7, False)],
        "hemi_engulph_homo",
    )


def test_merge_hemi_engulph_homo_by_one_on_right() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(0, 4, True)],
        [Deletion(0, 4, True), Deletion(4, 5, False)],
        "hemi_engulph_homo_by_one_on_right",
    )


def test_merge_hemi_engulph_homo_by_one_on_left() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(1, 5, True)],
        [Deletion(0, 1, False), Deletion(1, 5, True)],
        "hemi_engulph_homo_by_one_on_left",
    )


def test_merge_hemi_engulph_homo_by_one_on_both_sides() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(1, 4, True)],
        [Deletion(0, 1, False), Deletion(1, 4, True), Deletion(4, 5, False)],
        "hemi_engulph_homo_by_one_one_on_both_sides",
    )


def test_merge_homo_engulph_hemi_by_one_on_right() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(0, 4, False)],
        [Deletion(0, 5, True)],
        "homo_engulph_hemi_by_one_on_right",
    )


def test_merge_homo_engulph_hemi_by_one_on_left() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(1, 5, False)],
        [Deletion(0, 5, True)],
        "homo_engulph_hemi_by_one_on_left",
    )


def test_merge_homo_engulph_hemi_by_one_on_both_sides() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(1, 4, False)],
        [Deletion(0, 5, True)],
        "homo_engulph_hemi_by_one_on_both_sides",
    )


def test_merge_hemi_left_overlap_homo() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 7, False), Deletion(2, 9, True)],
        [Deletion(0, 2, False), Deletion(2, 9, True)],
        "hemi_left_overlap_homozygous",
    )


def test_merge_homo_engulf_hemi_same_start_pos() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(0, 7, True)],
        [Deletion(0, 7, True)],
        "homo_engulf_hemi_same_start_pos",
    )


def test_merge_hemi_engulf_homo_same_start_pos() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(0, 7, False)],
        [Deletion(0, 5, True), Deletion(5, 7, False)],
        "hemi_engulf_homo_same_start_pos",
    )


def test_merge_homo_engulf_hemi_same_end_pos() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(2, 7, False), Deletion(0, 7, True)],
        [Deletion(0, 7, True)],
        "homo_engulf_hemi_same_end_pos",
    )


def test_merge_hemi_engulf_homo_same_end_pos() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(2, 7, True), Deletion(0, 7, False)],
        [Deletion(0, 2, False), Deletion(2, 7, True)],
        "hemi_engulf_homo_same_end_pos",
    )


def test_merge_hemi_homo_same_size() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 7, True), Deletion(0, 7, False)],
        [Deletion(0, 7, True)],
        "hemi_homo_same_size",
    )


def test_single_homo_deletion() -> None:
    helper([Deletion(0, 7, True)], [Deletion(0, 7, True)], "single_homo_deletion")


def test_single_hemi_deletion() -> None:
    helper([Deletion(0, 7, False)], [Deletion(0, 7, False)], "single_hemi_deletion")


def test_multiple_non_overlapping_deletions() -> None:
    helper(
        [Deletion(0, 3, True), Deletion(5, 7, False), Deletion(10, 12, False)],
        [Deletion(0, 3, True), Deletion(5, 7, False), Deletion(10, 12, False)],
        "multiple_non_overlapping_deletions",
    )


def test_right_overlap_then_extra_deletion() -> None:
    helper(
        [Deletion(0, 5, True), Deletion(3, 7, False), Deletion(10, 12, False)],
        [Deletion(0, 5, True), Deletion(5, 7, False), Deletion(10, 12, False)],
        "right_overlap_then_extra_deletion",
    )


def test_two_consecutive_right_overlaps() -> None:
    helper(
        [
            Deletion(0, 5, True),
            Deletion(3, 7, False),
            Deletion(10, 15, False),
            Deletion(13, 17, True),
        ],
        [
            Deletion(0, 5, True),
            Deletion(5, 7, False),
            Deletion(10, 13, False),
            Deletion(13, 17, True),
        ],
        "two_consecutive_right_overlaps",
    )


def test_hemi_engulf_homo_split_followed_by_single_deletion() -> None:
    helper(
        [Deletion(0, 7, False), Deletion(2, 4, True), (Deletion(10, 12, True))],
        [
            Deletion(0, 2, False),
            Deletion(2, 4, True),
            Deletion(4, 7, False),
            Deletion(10, 12, True),
        ],
        "hemi_engulf_homo_split_followed_by_single_deletion",
    )


def test_homo_engulf_hemi_followed_by_single_deletion() -> None:
    helper(
        [Deletion(0, 7, True), Deletion(2, 4, False), (Deletion(10, 12, True))],
        [Deletion(0, 7, True), Deletion(10, 12, True)],
        "homo_engulf_hemi_followed_by_single_deletion",
    )


# Same zygosity combinations
def test_separate_homo_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(6, 11, True)],
        [Deletion(0, 5, True), Deletion(6, 11, True)],
        "separate_homo_deletions",
    )


def test_separate_hemi_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(6, 11, False)],
        [Deletion(0, 5, False), Deletion(6, 11, False)],
        "separate_hemi_deletions",
    )


def test_consecutive_homo_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 1, True), Deletion(1, 2, True)],
        [(Deletion(0, 2, True))],
        "consecutive_homo_deletions",
    )


def test_consecutive_hemi_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 1, False), Deletion(1, 2, False)],
        [(Deletion(0, 2, False))],
        "consecutive_hemi_deletions",
    )


def test_duplicate_homo_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [
            Deletion(0, 1, True),
            Deletion(0, 1, True),
            Deletion(0, 1, True),
            Deletion(0, 1, True),
        ],
        [Deletion(0, 1, True)],
        "duplicate_homo_deletions",
    )


def test_duplicate_hemi_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [
            Deletion(0, 1, False),
            Deletion(0, 1, False),
            Deletion(0, 1, False),
            Deletion(0, 1, False),
        ],
        [Deletion(0, 1, False)],
        "duplicate_hemi_deletions",
    )


def test_hemi_adjacent_homo() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(5, 10, True)],
        [Deletion(0, 5, False), Deletion(5, 10, True)],
        "hemi_adjacent_homo",
    )


def test_homo_adjacent_hemi() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(5, 10, True)],
        [Deletion(0, 5, False), Deletion(5, 10, True)],
        "homo_adjacent_hemi",
    )


def test_left_overlap_homo_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(3, 8, True)],
        [Deletion(0, 8, True)],
        "left_overlap_homo_deletions",
    )


def test_right_overlap_homo_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(3, 8, True), Deletion(0, 5, True)],
        [Deletion(0, 8, True)],
        "right_overlap_homo_deletions",
    )


def test_left_overlap_hemi_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(3, 8, False)],
        [Deletion(0, 8, False)],
        "left_overlap_hemi_deletions",
    )


def test_right_overlap_hemi_deletions() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(3, 8, False), Deletion(0, 5, False)],
        [Deletion(0, 8, False)],
        "right_overlap_hemi_deletions",
    )


def test_hemi_engulph_hemi_deletion() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, False), Deletion(2, 5, False)],
        [Deletion(0, 10, False)],
        "hemi_engulph_hemi_deletion",
    )


def test_homo_engulph_homo_deletion() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, True), Deletion(2, 5, True)],
        [Deletion(0, 10, True)],
        "homo_engulph_homo_deletion",
    )


def test_same_start_different_end_homo() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, True), Deletion(0, 8, True)],
        [Deletion(0, 8, True)],
        "same_start_different_end_homo",
    )


def test_same_end_different_start_homo() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(2, 8, True), Deletion(0, 8, True)],
        [Deletion(0, 8, True)],
        "same_end_different_start_homo",
    )


def test_same_start_different_end_hemi() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 5, False), Deletion(0, 8, False)],
        [Deletion(0, 8, False)],
        "same_start_different_end_hemi",
    )


def test_same_end_different_start_hemi() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(2, 8, False), Deletion(0, 8, False)],
        [Deletion(0, 8, False)],
        "same_end_different_start_hemi",
    )


def test_hemi_engulf_two_homos() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, False), Deletion(2, 4, True), Deletion(6, 8, True)],
        [
            Deletion(0, 2, False),
            Deletion(2, 4, True),
            Deletion(4, 6, False),
            Deletion(6, 8, True),
            Deletion(8, 10, False),
        ],
        "hemi_engulf_two_homos",
    )


def test_hemi_engulf_two_homos_one_at_start() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, False), Deletion(0, 3, True), Deletion(6, 8, True)],
        [
            Deletion(0, 3, True),
            Deletion(3, 6, False),
            Deletion(6, 8, True),
            Deletion(8, 10, False),
        ],
        "hemi_engulf_two_homos_one_at_start",
    )


def test_hemi_engulf_two_homos_one_at_end() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, False), Deletion(2, 4, True), Deletion(7, 10, True)],
        [
            Deletion(0, 2, False),
            Deletion(2, 4, True),
            Deletion(4, 7, False),
            Deletion(7, 10, True),
        ],
        "hemi_engulf_two_homos_one_at_end",
    )


def test_hemi_engulf_two_homos_at_both_ends() -> None:
    both_forward_and_reverse_merge_helper(
        [Deletion(0, 10, False), Deletion(0, 3, True), Deletion(7, 10, True)],
        [
            Deletion(0, 3, True),
            Deletion(3, 7, False),
            Deletion(7, 10, True),
        ],
        "hemi_engulf_two_homos_at_both_ends",
    )


def test_hemi_engulf_three_homos() -> None:
    helper(
        [
            Deletion(0, 20, False),
            Deletion(2, 5, True),
            Deletion(8, 11, True),
            Deletion(14, 17, True),
        ],
        [
            Deletion(0, 2, False),
            Deletion(2, 5, True),
            Deletion(5, 8, False),
            Deletion(8, 11, True),
            Deletion(11, 14, False),
            Deletion(14, 17, True),
            Deletion(17, 20, False),
        ],
        "hemi_engulf_three_homos",
    )


def test_three_overlapping_homos() -> None:
    helper(
        [Deletion(0, 5, True), Deletion(3, 8, True), Deletion(6, 11, True)],
        [Deletion(0, 11, True)],
        "three_overlapping_homos",
    )


def test_three_overlapping_hemis() -> None:
    helper(
        [Deletion(0, 5, False), Deletion(3, 8, False), Deletion(6, 11, False)],
        [Deletion(0, 11, False)],
        "three_overlapping_hemis",
    )


def test_hemi_homo_hemi_chain() -> None:
    helper(
        [Deletion(0, 5, False), Deletion(2, 8, True), Deletion(6, 11, False)],
        [Deletion(0, 2, False), Deletion(2, 8, True), Deletion(8, 11, False)],
        "hemi_homo_hemi_chain",
    )


def test_homo_hemi_homo_chain() -> None:
    helper(
        [Deletion(0, 5, True), Deletion(2, 8, False), Deletion(6, 11, True)],
        [Deletion(0, 5, True), Deletion(5, 6, False), Deletion(6, 11, True)],
        "homo_hemi_homo_chain",
    )


def test_chain_where_merge_changes_third_overlap() -> None:
    helper(
        [Deletion(0, 6, True), Deletion(3, 12, False), Deletion(9, 15, True)],
        [Deletion(0, 6, True), Deletion(6, 9, False), Deletion(9, 15, True)],
        "chain_where_merge_changes_third_overlap",
    )


def test_two_homos_spanning_hemi() -> None:
    helper(
        [Deletion(0, 5, True), Deletion(4, 10, True), Deletion(2, 7, False)],
        [Deletion(0, 10, True)],
        "two_homos_spanning_hemi",
    )


def test_duplicate_homo_and_hemi_same_position() -> None:
    helper(
        [
            Deletion(0, 5, True),
            Deletion(0, 5, True),
            Deletion(0, 5, False),
            Deletion(0, 5, False),
        ],
        [Deletion(0, 5, True)],
        "duplicate_homo_and_hemi_same_position",
    )


def test_two_homos_with_hemi_interleaved() -> None:
    helper(
        [Deletion(0, 10, True), Deletion(3, 6, False), Deletion(5, 15, True)],
        [Deletion(0, 15, True)],
        "two_homos_with_hemi_interleaved",
    )


def test_complex_alternating_homo_hemi_cascade() -> None:
    helper(
        [
            Deletion(0, 10, True),
            Deletion(5, 15, False),
            Deletion(12, 22, True),
            Deletion(18, 28, False),
            Deletion(25, 35, True),
        ],
        [
            Deletion(0, 10, True),
            Deletion(10, 12, False),
            Deletion(12, 22, True),
            Deletion(22, 25, False),
            Deletion(25, 35, True),
        ],
        "complex_alternating_homo_hemi_cascade",
    )


def test_complex_large_hemi_with_many_homos_inside() -> None:
    helper(
        [
            Deletion(0, 50, False),
            Deletion(5, 10, True),
            Deletion(15, 20, True),
            Deletion(25, 30, True),
            Deletion(35, 40, True),
        ],
        [
            Deletion(0, 5, False),
            Deletion(5, 10, True),
            Deletion(10, 15, False),
            Deletion(15, 20, True),
            Deletion(20, 25, False),
            Deletion(25, 30, True),
            Deletion(30, 35, False),
            Deletion(35, 40, True),
            Deletion(40, 50, False),
        ],
        "complex_large_hemi_with_many_homos_inside",
    )


def test_alternating_hemi_homo_pyramid() -> None:
    helper(
        [
            Deletion(15, 85, True),
            Deletion(20, 80, False),
            Deletion(25, 75, True),
            Deletion(30, 70, False),
            Deletion(35, 65, True),
            Deletion(40, 60, False),
            Deletion(45, 55, True),
        ],
        [
            Deletion(15, 85, True),
        ],
        "alternating_hemi_homo_pyramid",
    )


def test_empty() -> None:
    helper([], [], "empty")


def test_with_backwards_homozygous_intervals() -> None:
    with pytest.raises(InvertedDeletionError) as excinfo:
        merge_deletion_set(DeletionSet([Deletion(30, 15, True)], "test", "chr1"))
    assert (
        str(excinfo.value)
        == "Deletion set for 'test' in chromosome 'chr1' contains inverted homozygous deletion of range 30 - 15"
    )


def test_with_backwards_hemizygous_intervals() -> None:
    with pytest.raises(InvertedDeletionError) as excinfo:
        merge_deletion_set(DeletionSet([Deletion(30, 15, False)], "test", "chr1"))
    assert (
        str(excinfo.value)
        == "Deletion set for 'test' in chromosome 'chr1' contains inverted hemizygous deletion of range 30 - 15"
    )

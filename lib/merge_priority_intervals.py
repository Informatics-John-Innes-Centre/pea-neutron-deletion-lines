"""Generic algorithm for merging a set of intervals based on their priority"""

from dataclasses import dataclass
from collections import Counter
from heapq import heappush_max, heappop_max
from lib.deletion import Deletion, DeletionSet, InvertedDeletionError


@dataclass
class Interval:
    start: int
    end: int
    priority: int


class InvertedIntervalError(ValueError):
    def __init__(self, interval: Interval):
        self.interval = interval
        super().__init__(f"Inverted interval: [{interval.start}, {interval.end}]")


@dataclass
class Event:
    pos: int
    priority: int
    start: bool


def get_events_from_intervals(intervals: list[Interval]) -> list[Event]:
    events: list[Event] = []
    for interval in intervals:
        if interval.start > interval.end:
            raise InvertedIntervalError(interval)
        events.append(Event(interval.start, interval.priority, True))
        events.append(Event(interval.end, interval.priority, False))
    events.sort(key=lambda event: (event.pos, not event.start))
    return events


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    # Merge intervals with priority
    events = get_events_from_intervals(intervals)
    active: Counter[int] = Counter()
    result: list[Interval] = []
    heap: list[int] = []
    current_start = None
    previous_highest_priority = None
    for event in events:
        if event.start:
            # On start, pushing priority onto heap and keep counter
            active[event.priority] += 1
            heappush_max(heap, event.priority)
        else:
            if active[event.priority] > 1:
                active[event.priority] -= 1
            else:
                del active[event.priority]
            # Lazy delete: clean up heap top for priorities that are no longer active
            while heap and heap[0] not in active:
                heappop_max(heap)
        # heap[0] is now the highest priority active interval (most negative = highest)
        current_highest_priority = heap[0] if heap else None
        if current_highest_priority != previous_highest_priority:
            if previous_highest_priority is not None and current_start != event.pos:
                if current_start is None:
                    raise ValueError("Current start should not be None here.")
                result.append(
                    Interval(current_start, event.pos, previous_highest_priority)
                )
            if current_highest_priority is not None:
                current_start = event.pos
            previous_highest_priority = current_highest_priority
    return result


def merge_deletion_set(deletion_set: DeletionSet) -> DeletionSet:
    intervals = [
        Interval(deletion.start, deletion.end, 1 if deletion.homozygous else 0)
        for deletion in deletion_set.deletions
    ]

    try:
        merged_intervals = merge_intervals(intervals)
    except InvertedIntervalError as e:
        raise InvertedDeletionError(
            Deletion(e.interval.start, e.interval.end, e.interval.priority == 1),
            deletion_set.chromosome,
            deletion_set.line_sub,
        ) from e

    return DeletionSet(
        [
            Deletion(
                interval.start, interval.end, True if interval.priority == 1 else False
            )
            for interval in merged_intervals
        ],
        deletion_set.line_sub,
        deletion_set.chromosome,
    )

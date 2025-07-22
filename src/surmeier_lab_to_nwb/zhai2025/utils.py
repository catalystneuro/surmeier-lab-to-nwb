"""Utility functions for NWB conversion."""

from pynwb.base import TimeSeries


def time_series_duration(time_series: TimeSeries) -> float:
    """
    Calculate the duration of a TimeSeries using its timing information.

    If the TimeSeries has explicit timestamps, uses those to calculate duration.
    If it has starting_time and rate, calculates duration from data shape and rate.

    Parameters
    ----------
    time_series : TimeSeries
        The TimeSeries object to calculate duration for

    Returns
    -------
    float
        Duration in seconds
    """
    if hasattr(time_series, "timestamps") and time_series.timestamps is not None:
        # Use explicit timestamps
        timestamps = time_series.timestamps[:]
        return float(timestamps[-1] - timestamps[0])
    elif hasattr(time_series, "starting_time") and hasattr(time_series, "rate"):
        # Use starting_time and rate
        data_length = len(time_series.data)
        return (data_length - 1) / time_series.rate
    else:
        raise ValueError("TimeSeries must have either timestamps or starting_time+rate")

def is_gc_in_range(sequence: str, gc_bounds: int or float or tuple(int or float)) -> bool:
    """
    Determines if the GC content of a DNA sequence falls within a specified range.

    Parameters:
        seq (str): A DNA sequence as a string.
        gc_range (Union[float, Tuple[float, float]]): A single value or tuple of two values representing the desired
            range of GC content. If a single value is given, the function will check if the GC content is equal to that
            value. If a tuple is given, the function will check if the GC content is within the inclusive range of the
            two values.

    Returns:
        bool: True if the GC content of the sequence falls within the specified range, False otherwise.
    """
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    gc_count = (
        sequence.count("G")
        + sequence.count("C")
        + sequence.count("g")
        + sequence.count("c")
    )
    gc_percent = gc_count * 100 / len(sequence)
    return gc_bounds[0] <= gc_percent <= gc_bounds[1]


def is_length_in_range(sequence: str, length_bounds: int) -> bool:
    """
    Determines if the length of a DNA sequence falls within a specified range.

    Parameters:
        seq (str): A DNA sequence as a string.
        length_range (Union[int, Tuple[int, int]]): A single value or tuple of two values representing the desired
            range of sequence length. If a single value is given, the function will check if the length of the sequence
            is equal to that value. If a tuple is given, the function will check if the length of the sequence is within
            the inclusive range of the two values.

    Returns:
        bool: True if the length of the sequence falls within the specified range, False otherwise.
    """
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    return length_bounds[0] <= len(sequence) <= length_bounds[1]


def is_quality_in_range(quality: str, quality_threshold: int or float) -> bool:
    """
    Determines if the average quality score of a DNA sequence meets or exceeds a specified threshold.

    Parameters:
        quality_str (str): A string representing the quality scores of a DNA sequence.
        quality_threshold (Union[int, float]): An integer or float representing the desired quality threshold.

    Returns:
        bool: True if the average quality score of the sequence meets or exceeds the specified threshold, False otherwise.
    """
    quality_scores = []
    for symbol in quality:
        quality_scores.append(ord(symbol) - 33)
    average_quality = sum(quality_scores) / len(quality)
    return average_quality >= quality_threshold


def fastq_filter(
    seqs: dict{str, tuple(str) or list(str)},
    gc_bounds: int or float or tuple(int or float) or list[int or float] = (0, 100),
    length_bounds: int or tuple[int] or list[int] = (0, 2**32),
    quality_threshold: int or float = 0,
) -> dict:
    """
    Filters a dictionary of DNA sequences and their quality scores based on specified criteria.

    Parameters:
        seq_dict (Dict[str, Tuple[str, str]]): A dictionary mapping sequence names to tuples containing the DNA sequence
            as a string and its quality scores as a string.
        gc_range (Optional[Union[float, Tuple[float, float]]]): A single value or tuple of two values representing
            the desired range of GC content. If a single value is given, the function will check if the GC content is
            equal to that value. If a tuple is given, the function will check if the GC content is within the inclusive
            range of the two values.
        length_range (Optional[Union[int, Tuple[int, int]]]): A single value or tuple of two values representing the
            desired range of sequence length. If a single value is given, the function will check if the length of the
            sequence is equal to that value. If a tuple is given, the function will check if the length of the sequence
            is within the inclusive range of the two values.
        quality_threshold (Optional[Union[int, float]]): An integer or float representing the desired quality threshold.

    Returns:
        Dict[str, Tuple[str, str]]: A dictionary containing only the DNA sequences and their quality scores that pass
            all specified filters.
    """
    fastq_dictionary = {}
    for sequence_name in seqs.keys():
        sequence = seqs[sequence_name][0]
        quality = seqs[sequence_name][1]
    if is_gc_in_range(sequence=sequence, gc_bounds=gc_bounds):
        if is_length_in_range(sequence, length_bounds):
            if is_quality_in_range(quality, quality_threshold):
                fastq_dictionary[sequence_name] = seqs[sequence_name]
    return fastq_dictionary

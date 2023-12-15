from bioinf_modules import dna_rna_tools
from bioinf_modules import fastq_filter
from bioinf_modules import protein_tools


def run_dna_rna_tools(*args):
    procedure, sequences = args[-1], args[:-1]
    if procedure not in dna_rna_tools.PROCEDURES_DICT.keys():
        raise ValueError("Wrong procedure")
    rna_dna_tools_result_list = []
    for sequence in sequences:
        dna_rna_tools.check_user_input(sequence)
        rna_dna_tools_result_list.append(
            dna_rna_tools.PROCEDURES_DICT[procedure](sequence)
        )
    if len(rna_dna_tools_result_list) == 1:
        return rna_dna_tools_result_list[0]
    else:
        return rna_dna_tools_result_list


def run_protein_tools(sequences: (str, tuple[str] or list[str]), **kwargs: str) -> dict:
    """
    Process protein sequence by one of the developed tools.\n
    Run one procedure at a time:
    - Search for conserved amino acids residues in protein sequence
    - Search for alternative frames in a protein sequences
    - Convert protein sequences to RNA or DNA sequences
    - Reverse the protein sequences from one-letter to three-letter format and vice-versa
    - Define molecular weight of the protein sequences

    All functions except *search_for_alt_frames* are letter case sensitive\n
    If only one sequence provided - *sequences* can be string.\n
    If more - please provide *sequences* as list or tuple.\n
    Provide protein sequence in one letter code.\n
    You can obtain one letter code from three letter code with *three_one_letter_code*\n
    If more information needed please see README or desired docstring

    Arguments:
    - sequences (str, list[str] or tuple[str]): sequences to process
    - procedure (str]: desired procedure:
        - "search_for_motifs"
        - "search_for_alt_frames"
        - "convert_to_nucl_acids"
        - "three_one_letter_code"
        - "define_molecular_weight"

    For "search_for_motif" procedure provide:
    - motif (str): desired motif to check presense in every given sequence\n
            Example: motif = "GA"
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
            Example: overlapping = False

    For "search_for_alt_frames" procedure provide:
    - alt_start_aa (str): the name of an amino acid that is encoded by alternative start codon (Optional)\n
            Example: alt_start_aa = "I"

    For "convert_to_nucl_acids" procedure provide:
    - nucl_acids (str): the nucleic acid to convert to\n
            Example: nucl_acids = "RNA"\n
                           nucl_acids = "DNA"\n
                           nucl_acids = "both"

    Return:
    - dict: Dictionary with processed sequences. Depends on desired tool\n
            Please see Readme or desired docstring
    """
    procedure_arguments, procedure = protein_tools.check_and_parse_user_input(
        sequences, **kwargs
    )
    return protein_tools.PROTEINS_PROCEDURES_TO_FUNCTIONS[procedure](
        **procedure_arguments
    )


def run_fastq_filter(
    seqs: dict[str, tuple[str] | list[str]],
    gc_bounds: int | float | tuple[int | float] | list[int | float] = (0, 100),
    length_bounds: int | tuple[int] | list[int] = (0, 2**32),
    quality_threshold: int | float = 0,
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
    if fastq_filter.is_gc_in_range(sequence=sequence, gc_bounds=gc_bounds):
        if fastq_filter.is_length_in_range(sequence, length_bounds):
            if fastq_filter.is_quality_in_range(quality, quality_threshold):
                fastq_dictionary[sequence_name] = seqs[sequence_name]
    return fastq_dictionary

run_protein_tools(['mNYQTMSPYYDMId'], procedure='search_for_alt_frames')  # {'mNYQTMSPYYDMId': ['MSPYYDMId']}


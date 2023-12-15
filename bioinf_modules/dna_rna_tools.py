def check_user_input(sequence):
    if not all(i in "".join(COMPLEMENT_DICT.keys()) for i in sequence):
        raise ValueError("Invalid sequence given")
    if "T" in sequence.upper() and "U" in sequence.upper():
        raise ValueError("DNA-RNA mix sequence given")


def reverse(sequence):
    return sequence[::-1]


def is_dna(sequence):
    return "T" in sequence or "t" in sequence


def transcribe(sequence):
    if is_dna:
        return sequence.replace("T", "U").replace("t", "u")
    else:
        return sequence.replace("U", "T").replace("u", "t")


COMPLEMENT_DICT = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "T": "A",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
    "t": "a",
}


def complement(sequence):
    if is_dna:
        COMPLEMENT_DICT["A"] = "T"
        COMPLEMENT_DICT["a"] = "t"
    complement_seq = ""
    for letter in sequence:
        complement_seq = complement_seq + COMPLEMENT_DICT[letter]
    return complement_seq


def reverse_complement(sequence):
    if is_dna:
        COMPLEMENT_DICT["A"] = "T"
        COMPLEMENT_DICT["a"] = "t"


PROCEDURES_DICT = {
    "transcribe": transcribe,
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
}

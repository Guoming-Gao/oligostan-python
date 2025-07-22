# filters.py
try:
    from Bio.SeqUtils import gc_fraction
except ImportError:
    from Bio.SeqUtils import GC

    def gc_fraction(seq):
        return GC(seq) / 100


def is_ok_4_pnas_filter(seq, filter_to_be_used=[1, 2, 3, 4, 5]):
    """Exact translation of isOk4PNASFilter from R"""
    results = []

    if 1 in filter_to_be_used:
        results.append(is_it_ok_4_a_comp(seq))
    if 2 in filter_to_be_used:
        results.append(is_it_ok_4_a_stack(seq))
    if 3 in filter_to_be_used:
        results.append(is_it_ok_4_c_comp(seq))
    if 4 in filter_to_be_used:
        results.append(is_it_ok_4_c_stack(seq))
    if 5 in filter_to_be_used:
        results.append(is_it_ok_4_c_spec_stack(seq))

    return all(results)


def is_it_ok_4_a_comp(seq):
    """PNAS Rule 1: Adenine content < 28%"""
    a_count = seq.upper().count("A")
    return (a_count / len(seq)) < 0.28


def is_it_ok_4_a_stack(seq):
    """PNAS Rule 2: No AAAA runs"""
    return "AAAA" not in seq.upper()


def is_it_ok_4_c_comp(seq):
    """PNAS Rule 3: Cytosine content between 22-28%"""
    c_count = seq.upper().count("C")
    c_comp = c_count / len(seq)
    return 0.22 < c_comp < 0.28


def is_it_ok_4_c_stack(seq):
    """PNAS Rule 4: No CCCC runs"""
    return "CCCC" not in seq.upper()


def is_it_ok_4_c_spec_stack(seq):
    """PNAS Rule 5: No 6-nt windows with >50% cytosine"""
    seq = seq.upper()
    for i in range(len(seq) - 5):
        window = seq[i : i + 6]
        c_percent = window.count("C") / 6
        if c_percent > 0.5:
            return False
    return True


def is_ok_4_gc_filter(seq, min_gc=0.4, max_gc=0.6):
    """GC content filter matching R script logic"""
    gc_content = gc_fraction(seq.upper())
    return min_gc <= gc_content <= max_gc

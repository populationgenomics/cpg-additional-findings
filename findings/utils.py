"""
utility methods and constants
"""


import json
from typing import Any

from cpg_utils import to_path


AMINO_ACIDS = {
    'C': 'Cys',
    'D': 'Asp',
    'S': 'Ser',
    'Q': 'Gln',
    'K': 'Lys',
    'I': 'Ile',
    'P': 'Pro',
    'T': 'Thr',
    'F': 'Phe',
    'N': 'Asn',
    'G': 'Gly',
    'H': 'His',
    'L': 'Leu',
    'R': 'Arg',
    'W': 'Trp',
    'A': 'Ala',
    'V': 'Val',
    'E': 'Glu',
    'Y': 'Tyr',
    'M': 'Met',
}


def read_json_from_path(bucket_path: str) -> Any:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with to_path(bucket_path).open() as handle:
        return json.load(handle)

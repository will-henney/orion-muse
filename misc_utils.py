"""
Miscellaneous routines
"""
import sys

py3 = sys.version_info >= (3, 0)

if py3:
    TRANSLATION_TABLE = {
        ord(" "): '_',
        ord(","): None,
        ord("["): None,
        ord("]"): None,
        ord("("): '_',
        ord(")"): None,
    }
else:
    import string
    TRANSLATION_TABLE = string.maketrans(' (', '__')

def sanitize_string(s):
    """Make string suitable to be a unix filename
Remove and/or substitute problematic characters via TRANSLATION_TABLE"""
    if py3:
        return s.translate(TRANSLATION_TABLE)
    else:
        return s.translate(TRANSLATION_TABLE, ',[])')

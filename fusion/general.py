'''=========================================================

	F U N C T I O N S

========================================================='''

'''---------------------------------------------------------
#
#   Strip a series of characters from a string
#
#   Input:    A string to strip, A list of characters to remove
#   Output:   The stripped string.
#
---------------------------------------------------------'''

def stripStringChars(to_strip, to_remove):

    for remove_this in to_remove:
        to_strip = to_strip.strip(remove_this)

    return to_strip

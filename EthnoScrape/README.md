EthnoScrape
====================

# Purpose

This is a very simple script to "scrape" the Ethnologue website (using
BeautifulSoup) to get the genetic classification for many of the languages in
WALS.  The ISO codes are used to map between languages in the two databases.
This is necessary because, unlike WALS, Ethnologue do not, as far as I can
tell, make all of their data available to download in a convenient fashion.

# Disclaimer

I'm not sure how the Ethnologue folks would necessarily feel about this, but
the behaviour of the script is not in violation of Ethnologue's robots.txt
(as of early March 2013).  Their robot file is moderately detailed, so this
seems unlikely to be an oversight.

# Files

`iso_codes.txt` is a plain text file.  Each line contains one or more language
ISO codes.  Lines with multiple codes represent (I think) multiple codes for
the same language (which seems silly).  Multiple codes on one line ar
separated by spaces.

`ethnoscrape.py` is a Python script which takes no arguments.  It expects
`iso_codes.txt` to exist and to be as described above.  It produces as output
the file `classifications.txt` described below.  The output is relatively
verbose, so you can track the progress.  It waits 10 seconds between
consecutive hits to the Ethnologue website so as not to cause them undue server
load.  There are six hundred and something ISO codes in `iso_codes.txt`, so the
whole scrape takes over an hour.

`classifications.txt` is the output of running `ethnoscrape.py`.  It contains
one line for every ISO code which was successfully looked up.  Each line is
the iso code, followed by a colon, a space and then the classification string,
which is a comma-separated list of nested language groups.

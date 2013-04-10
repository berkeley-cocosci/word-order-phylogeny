WALS2SQL
====================

# Purpose

This script takes the WALS data, in the standard directory of .csv files format
that WALS makes available for download, and turns it into an SQLite database,
to make it easy access WALS data from within other Python code.

# Files

`wals_data` is a directory of .csv files containing all of the data in WALS.

`wals2sql.py` is a Python script which takes no arguments.  It expects the
`wals_data` directory to exist and to be as described above.  It produces as
output the file `wals.db` as described below.

`wals.db` is the output of running `wals2sql.py`.  It is an SQLite database
file containing all of the data in WALS.

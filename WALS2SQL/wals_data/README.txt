
WALS Data
=========

WALS Online data export created 2011-10-23 23:27:51.979916 from http://wals.info/

Contents
--------

File Name                                              Size
wals_data/datapoints.csv                             606377
wals_data/values.csv                                  67660
wals_data/features.csv                                 7312
wals_data/languages.csv                              169264

File encoding: utf8


Data Structure
--------------

datapoints.csv contains a matrix of all value assignments for features in WALS.

Rows give the value assignments for a particular language identified by its
WALS code given in the first column,

Columns give the value assignments for a particular feature identified by its
numeric identifier given in the first row.

Descriptions of the values can be looked up in values.csv, using the numeric
feature id and the numeric value id given in the datapoints matrix.

Additional data for languages can be looked up using the WALS code in
languages.csv

Additional data for features can be looked up using the numeric feature id in
features.csv

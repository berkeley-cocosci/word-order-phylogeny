#!/usr/bin/env python
from codecs import open
import csv
import sqlite3

#from make_dist_matrix import *

class Language:

    def __init__(self):
	self.code = "FOO"
	self.name = "BAR"
	self.data = {}
		
def language_from_wals_code_factory(code, conn, cursor):
    lang = Language()
    cursor.execute('''SELECT * FROM dense_languages WHERE wals_code=?''',(code,))
    results = cursor.fetchone()
    lang.code = results[0]
    lang.name = results[1].replace(" ","_").replace("(","").replace(")","")
    lang.data = {}
    lang.data["location"] = (float(results[2]),float(results[3]))
    lang.data["genus"] = results[4]
    lang.data["family"] = results[5]
    lang.data["subfamily"] = results[6]
    cursor.execute('''SELECT features.name, dense_data.value_id FROM dense_data INNER JOIN features on dense_data.feature_id = features.id WHERE wals_code=?''',(code,))
    for x in cursor.fetchall():
        name, value = x
        lang.data[name] = value
    return lang

def create_languages_table(conn, cursor):
    cursor.execute('''CREATE TABLE languages(
                wals_code TEXT PRIMARY KEY,
                name TEXT,
                latitude REAL,
                longitude REAL,
                genus TEXT,
                family TEXT,
                subfamily TEXT,
                iso_codes TEXT)''')

def create_features_table(conn, cursor):
    cursor.execute('''CREATE TABLE features(
                id TEXT PRIMARY KEY,
                name TEXT)''')

def create_values_table(conn, cursor):
    cursor.execute('''CREATE TABLE valuess(
                feature_id TEXT,
                value_id INTEGER,
                short_desc TEXT,
                long_desc TEXT)''')

def create_data_table(conn, cursor):
    cursor.execute('''CREATE TABLE data_points(
                wals_code TEXT,
                feature_id TEXT,
                value_id INTEGER,
                FOREIGN KEY(wals_code) REFERENCES languages(wals_code),
                FOREIGN KEY(feature_id) REFERENCES features(id))''')
#                FOREIGN KEY(value_id) REFERENCES valuess(value_id))''')

def create_tables(conn, cursor):
    create_languages_table(conn, cursor)
    create_features_table(conn, cursor)
    create_values_table(conn, cursor)
    create_data_table(conn, cursor)

def parse_wals_data(wals_dir):
    # Build a data structure representing all the WALS data
    # This structure is a dictionary with one key/value pair per file of WALS data
    # The dictionary keys are the filenames (e.g. "languages", "features", etc.)
    # The dictionary values are lists of dictionaries
    # Each dictionary in a main dictionary value list represents one line from the file
    # The keys are the column headers, the values are the column contents
    data = {}
    filenames = "languages features values datapoints".split()
    for filename in filenames:
        x = []
        fp = open(wals_dir+"/"+filename+".csv", "r")
        reader = csv.DictReader(fp)
        for line in reader:
            x.append(line)
        fp.close()
        data[filename] = x
    return data

def populate_languages_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO languages VALUES (?, ?, ?, ?, ?, ?, ?,?)""", (unicode(datum["wals code"],"utf8"), unicode(datum["name"],"utf8"), float(datum["latitude"]), float(datum["longitude"]), unicode(datum["genus"],"utf8"), unicode(datum["family"],"utf8"), unicode(datum["subfamily"],"utf8"), unicode(datum["iso codes"],"utf8")))
def populate_features_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO features VALUES (?, ?)""", (datum["id"], unicode(datum["name"],"utf8")))

def populate_values_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO valuess VALUES (?, ?, ?, ?)""", (datum["feature_id"], int(datum["value_id"]), unicode(datum["description"],"utf8"), unicode(datum["long description"],"utf8")))

def populate_data_table(data, conn, cursor):
    for datum in data:
        keys = datum.keys()
        keys.remove("wals_code")
        for key in keys:
            wals_code = unicode(datum["wals_code"],"utf8")
            feature_id = key
            if datum[key] == "":
                value_id = None
            else:
                value_id = datum[key]
            cursor.execute("""INSERT INTO data_points VALUES (?, ?, ?)""", (wals_code, feature_id, value_id))

def insert_data(wals_dir, conn, cursor):
    wals_data = parse_wals_data(wals_dir)
    populate_languages_table(wals_data["languages"],conn, cursor)
    conn.commit()
    populate_features_table(wals_data["features"],conn, cursor)
    conn.commit()
    populate_values_table(wals_data["values"],conn, cursor)
    conn.commit()
    populate_data_table(wals_data["datapoints"],conn, cursor)    
    conn.commit()

def make_count_tables(conn, cursor):
    cursor.execute('''CREATE TABLE langs_per_feature_counts(
            feature_id TEXT,
            count INTEGER,
            FOREIGN KEY(feature_id) REFERENCES features(id))''')
    
    cursor.execute('''CREATE TABLE features_per_lang_counts(
            wals_code TEXT,
            count INTEGER,
            FOREIGN KEY(wals_code) REFERENCES languages(wals_code))''')

def populate_count_tables(conn, cursor):
    cursor.execute('''SELECT id FROM features''')
    results = cursor.fetchall()
    for result in results:
        feature_id = result[0]
        cursor.execute('''SELECT COUNT(wals_code) FROM data_points WHERE feature_id=? AND value_id IS NOT NULL''',(feature_id,))
        count = int(cursor.fetchone()[0])
        cursor.execute('''INSERT INTO langs_per_feature_counts VALUES (?,?)''', (feature_id, count))

    cursor.execute('''SELECT wals_code FROM languages''')
    results = cursor.fetchall()
    for result in results:
        wals_code = result[0]
        cursor.execute('''SELECT COUNT(feature_id) FROM data_points WHERE wals_code=? AND value_id IS NOT NULL''',(wals_code,))
        count = int(cursor.fetchone()[0])
        cursor.execute('''INSERT INTO features_per_lang_counts VALUES (?,?)''', (wals_code, count))
        
def compute_counts(conn, cursor):
    make_count_tables(conn, cursor)
    populate_count_tables(conn, cursor)

def report(conn, cursor):
    cursor.execute('''SELECT name, count FROM features INNER JOIN langs_per_feature_counts ON features.id = langs_per_feature_counts.feature_id ORDER BY count DESC LIMIT 10''')
    for result in cursor.fetchall():
        print result

def create_dense_subset(conn, cursor):
    cursor.execute('''DROP TABLE dense_languages''')
    cursor.execute('''CREATE TABLE dense_languages(
                wals_code TEXT PRIMARY KEY,
                name TEXT,
                latitude REAL,
                longitude REAL,
                genus TEXT,
                family TEXT,
                subfamily TEXT,
                iso_codes TEXT)''')
    cursor.execute('''INSERT INTO dense_languages SELECT * FROM languages''')

    #cursor.execute('''DROP TABLE dense_features''')
    cursor.execute('''CREATE TABLE dense_features(
                id TEXT PRIMARY KEY,
                name TEXT)''')
  
    n = 17 
    # Get the n+1 densest features
    # Why n+1?  Because Basic Word Order is dense, but we don't want to use it for our purposes
    # So if we grab the n+1 densest and exclude BWO, we get the n densest *usable* features
    cursor.execute('''SELECT id FROM features INNER JOIN langs_per_feature_counts ON features.id = langs_per_feature_counts.feature_id ORDER BY count DESC LIMIT ?''', (n+1,))
    dense_features = cursor.fetchall()
    dense_features = [x[0] for x in dense_features]
    for fid in dense_features:
        cursor.execute('''INSERT INTO dense_features SELECT * FROM features WHERE id=?''',(fid,))
        cursor.execute('''SELECT wals_code FROM data_points WHERE feature_id=? and value_id IS NULL''', (fid,))
        langs_to_kill = cursor.fetchall()
        langs_to_kill = [x[0] for x in langs_to_kill]
        for lang in langs_to_kill:
            cursor.execute('''DELETE FROM dense_languages WHERE wals_code=?''', (lang,))
    cursor.execute('''SELECT COUNT(wals_code) FROM dense_languages''')
    dense_count = int(cursor.fetchone()[0])
    print "There are %d dense languages" % dense_count
    cursor.execute('''SELECT name, genus, family, subfamily FROM dense_languages''')
    for lang in cursor.fetchall():
        print lang

def instantiate_dense_language_objects(conn, cursor):
    languages = []
    cursor.execute('''DROP TABLE dense_data''')
    cursor.execute('''CREATE TABLE dense_data AS
        SELECT wals_code, feature_id, value_id FROM data_points
        WHERE wals_code IN (SELECT wals_code FROM dense_languages) AND feature_id IN (SELECT id FROM dense_features)''')
        
##    cursor.execute('''CREATE TABLE big_dense_table AS
##        SELECT dense_languages.wals_code, dense_features.id, data_points.value_id FROM
##        dense_languages CROSS JOIN dense_features CROSS JOIN data_points WHERE
##        dense_languages.wals_code = data_points.wals_code AND
##        dense_features.id = data_points.feature_id''') 
    cursor.execute('''SELECT wals_code FROM dense_languages''')
    dense_codes = [x[0] for x in cursor.fetchall()]
    i = 0
    for code in dense_codes:
        i += 1
        print "Working on language %d of %d..." % (i, len(dense_codes))
        languages.append(language_from_wals_code_factory(code, conn, cursor))
    return languages


def save_translate_file(languages):
    fp = open("translate_block", "w", "utf8")
    fp.write("translate\n");
    fp.write("%d %s" % (0, languages[0].name)) 
    for i in range(1,len(languages)):
        fp.write(",\n%d %s" % (i, languages[i].name)) 
    fp.write(";\n")
    fp.close()

def save_multistate_file(languages):
    fp = open("multistate_file", "w", "utf8")
    for i, lang in enumerate(languages):
	if len(lang.data.keys()) != 19:
		print "Skipping %s, because it's broken" % lang.name.encode("UTF-8")
		continue
	fp.write("%s	%d\n" % (lang.name, lang.data[u'Order of Subject, Object and Verb']))
    fp.close()


def main():
    conn = sqlite3.connect("walssql2.db")
    cursor = conn.cursor()

    # Enable foreign keys
    cursor.execute('''PRAGMA foreign_keys = ON;''')

    # Create tables
    print "Creating tables...",
    create_tables(conn, cursor)
    print "done."
    
    # Write data into tables
    print "Inserting WALS data into tables...",
    insert_data("wals_data", conn, cursor)
    print "done."
    
    # Compute statistics
    print "Computing feature and language counts",
    compute_counts(conn, cursor)
    print "done."
    
    # Save (commit) the changes
    conn.commit()

    create_dense_subset(conn, cursor)
    # Save (commit) the changes
    conn.commit()

#    languages = instantiate_dense_language_objects(conn, cursor)
    
    # We can also close the cursor if we are done with it
    cursor.close()
    conn.close()

#    save_translate_file(languages)
#    save_multistate_file(languages)

#    make_matrix_go_now(languages)

if __name__ == "__main__":
    main()

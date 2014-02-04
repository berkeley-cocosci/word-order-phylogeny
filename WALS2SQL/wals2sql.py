#!/usr/bin/env python
from codecs import open
import csv
import sqlite3
import sys
sys.path.append("../TreeBuilder/")
import generate_trees

bwo = "Order of Subject, Object and Verb"
#from make_dist_matrix import *

class Language:

    def __init__(self):
	self.code = "FOO"
	self.name = "BAR"
	self.data = {}

def create_languages_table(conn, cursor):
    cursor.execute('''CREATE TABLE IF NOT EXISTS languages(
                wals_code TEXT PRIMARY KEY,
                name TEXT,
                latitude REAL,
                longitude REAL,
                genus TEXT,
                family TEXT,
                subfamily TEXT,
                iso_codes TEXT)''')

def create_features_table(conn, cursor):
    cursor.execute('''CREATE TABLE IF NOT EXISTS features(
                id TEXT PRIMARY KEY,
                name TEXT)''')

def create_values_table(conn, cursor):
    cursor.execute('''CREATE TABLE IF NOT EXISTS values_(
                feature_id TEXT,
                value_id INTEGER,
                short_desc TEXT,
                long_desc TEXT)''')

def create_data_table(conn, cursor):
    cursor.execute('''CREATE TABLE IF NOT EXISTS data_points(
                wals_code TEXT,
                feature_id TEXT,
                value_id INTEGER,
                FOREIGN KEY(wals_code) REFERENCES languages(wals_code),
                FOREIGN KEY(feature_id) REFERENCES features(id))''')

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

def empty_tables(conn, cursor):
    cursor.execute("""DELETE FROM data_points""")
    cursor.execute("""DELETE FROM langs_per_feature_counts""")
    cursor.execute("""DELETE FROM features_per_lang_counts""")
    cursor.execute("""DELETE FROM languages""")
    cursor.execute("""DELETE FROM features""")
    cursor.execute("""DELETE FROM values_""")

def populate_languages_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO languages VALUES (?, ?, ?, ?, ?, ?, ?,?)""", (unicode(datum["wals code"],"utf8"), unicode(datum["name"],"utf8"), float(datum["latitude"]), float(datum["longitude"]), unicode(datum["genus"],"utf8"), unicode(datum["family"],"utf8"), unicode(datum["subfamily"],"utf8"), unicode(datum["iso codes"],"utf8")))
def populate_features_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO features VALUES (?, ?)""", (datum["id"], unicode(datum["name"],"utf8")))

def populate_values_table(data, conn, cursor):
    for datum in data:
        cursor.execute("""INSERT INTO values_ VALUES (?, ?, ?, ?)""", (datum["feature_id"], int(datum["value_id"]), unicode(datum["description"],"utf8"), unicode(datum["long description"],"utf8")))

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
    empty_tables(conn, cursor)
    populate_languages_table(wals_data["languages"],conn, cursor)
    populate_features_table(wals_data["features"],conn, cursor)
    populate_values_table(wals_data["values"],conn, cursor)
    populate_data_table(wals_data["datapoints"],conn, cursor)    

def make_count_tables(conn, cursor):
    cursor.execute('''CREATE TABLE IF NOT EXISTS langs_per_feature_counts(
            feature_id TEXT,
            count INTEGER,
            FOREIGN KEY(feature_id) REFERENCES features(id))''')

    cursor.execute('''CREATE TABLE IF NOT EXISTS features_per_lang_counts(
            wals_code TEXT,
            count INTEGER,
            FOREIGN KEY(wals_code) REFERENCES languages(wals_code))''')

    cursor.execute('''CREATE INDEX IF NOT EXISTS data_points_wals_code_index
            ON data_points(wals_code)''')
    cursor.execute('''CREATE INDEX IF NOT EXISTS data_points_feature_id_index
            ON data_points(feature_id)''')
    cursor.execute('''CREATE INDEX IF NOT EXISTS data_points_value_id_index
            ON data_points(value_id)''')

def populate_count_tables(conn, cursor):
    cursor.execute('''DELETE FROM langs_per_feature_counts''')
    cursor.execute('''SELECT id FROM features''')
    results = cursor.fetchall()
    for result in results:
        feature_id = result[0]
        cursor.execute('''SELECT COUNT(wals_code) FROM data_points WHERE feature_id=? AND value_id IS NOT NULL''',(feature_id,))
        count = int(cursor.fetchone()[0])
        cursor.execute('''INSERT INTO langs_per_feature_counts VALUES (?,?)''', (feature_id, count))

    cursor.execute('''DELETE FROM features_per_lang_counts''')
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

def compute_dense_features(conn, cursor, n=25):

    cursor.execute('''CREATE TABLE IF NOT EXISTS dense_features(
                id TEXT PRIMARY KEY,
                name TEXT)''')
    cursor.execute('''DELETE FROM dense_features''')
  
    # Get the n densest non-disqualified features
    # Excludes feaures which are either BWO or heavily correlated with BWO
    # BWO is 81A
    # We exclude 82A, 83A, 95A etc. because they are highly correlated
    # After excluding those, 81A is the densest feature.
    # So if we get the n+1 densest, we get the n denses valid features, plus BWO
    # TODO come up with a good reason why not to use BWO for tree building?
    cursor.execute('''SELECT id FROM features INNER JOIN langs_per_feature_counts ON features.id = langs_per_feature_counts.feature_id WHERE id NOT IN ("82A", "83A", "95A", "97A", "96A", "144A") ORDER BY count DESC LIMIT ?''', (n+1,))
    dense_features = cursor.fetchall()
    dense_features = [x[0] for x in dense_features]
    for fid in dense_features:
        cursor.execute('''INSERT INTO dense_features SELECT * FROM features WHERE id=?''', (fid,))

def get_dense_features(conn, cursor):
    cursor.execute('''SELECT name FROM dense_features''')
    return [x[0] for x in cursor.fetchall()]

def compute_dense_languages(conn, cursor):

    # Create or empty dense_languages table
    cursor.execute('''CREATE TABLE IF NOT EXISTS dense_languages(
                wals_code TEXT PRIMARY KEY,
                name TEXT,
                latitude REAL,
                longitude REAL,
                genus TEXT,
                family TEXT,
                subfamily TEXT,
                iso_codes TEXT)''')
    cursor.execute('''DELETE FROM dense_languages''')

    # Copy all languages into dense_languages
    cursor.execute('''INSERT INTO dense_languages SELECT * FROM languages''')

    # Now we delete from the dense_languages table those languages
    # which do not have data for all the dense features
    for fid in dense_features:
        cursor.execute('''DELETE FROM dense_languages WHERE wals_code IN
                (SELECT wals_code FROM data_points WHERE feature_id=? and value_id IS NULL)''', (fid,))

### Convenience functions

def get_languages_by_family(conn, cursor, family):
    ethnoclasses = generate_trees.load_ethnologue_classifications()
    cursor.execute("""SELECT wals_code FROM languages WHERE family=? AND iso_codes != ''""", (family,))
    codes = [code[0] for code in cursor.fetchall()]
    languages = map(lambda(x): language_from_wals_code(conn, cursor, x), codes)
    languages = map(lambda(x): generate_trees.apply_ethnoclass(x, ethnoclasses), languages)
    languages = filter(lambda(x): x.data["ethnoclass"].split(",")[0].strip() == x.data["family"], languages)
    languages = filter(lambda(x): x.data.get(bwo, None) not in (7,'7',None), languages)
    hierlengths = [len(x.data["ethnoclass"].split(",")[1:]) for x in languages ]
    return languages

def get_dense_languages_by_family(conn, cursor, family):
    ethnoclasses = generate_trees.load_ethnologue_classifications()
    cursor.execute("""SELECT wals_code FROM dense_languages WHERE family=? AND iso_codes != ''""", (family,))
    codes = [code[0] for code in cursor.fetchall()]
    languages = map(lambda(x): language_from_wals_code(conn, cursor, x), codes)
    languages = map(lambda(x): generate_trees.apply_ethnoclass(x, ethnoclasses), languages)
    languages = filter(lambda(x): x.data["ethnoclass"].split(",")[0].strip() == x.data["family"], languages)
    languages = filter(lambda(x): x.data.get(bwo, None) not in (7,'7',None), languages)
    hierlengths = [len(x.data["ethnoclass"].split(",")[1:]) for x in languages ]
    return languages

def language_from_wals_code(conn, cursor, code):
    lang = Language()
    cursor.execute('''SELECT * FROM languages WHERE wals_code=?''',(code,))
    results = cursor.fetchone()
    lang.code = results[0]
    lang.name = results[1].replace(" ","_").replace("(","").replace(")","")
    lang.data = {}
    lang.data["location"] = (float(results[2]),float(results[3]))
    lang.data["genus"] = results[4]
    lang.data["family"] = results[5]
    lang.data["subfamily"] = results[6]
    lang.data["iso_codes"] = results[7]
    cursor.execute('''SELECT name, value_id FROM speedyfeatures WHERE wals_code=?''',(code,))
    for x in cursor.fetchall():
        name, value = x
        lang.data[name] = value
    return lang


def main():
    conn = sqlite3.connect("wals.db")
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
    
    conn.commit()

    # Find dense subset of language x feature matrix
    compute_dense_features(conn, cursor)
    compute_dense_languages(conn, cursor)
    
    conn.commit()
    cursor.close()
    conn.close()

if __name__ == "__main__":
    main()

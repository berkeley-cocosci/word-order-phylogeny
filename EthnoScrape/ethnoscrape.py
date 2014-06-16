#!/usr/bin/env python
import codecs
import time
import urllib

import BeautifulSoup as bs

def get_classification_from_iso(iso_code):
    url = "http://www.ethnologue.com/language/%s" % iso_code
    handle = urllib.urlopen(url)
    soup = bs.BeautifulSoup(handle)
    classidiv = soup.find("div", "field-name-language-classification-link")
    classi = classidiv.find("div", "field-item")
    return classi.text

def read_iso_codes_from_file(filename):
    fp = open(filename, "r")
    codes = []
    for line in fp:
        for code in line.split():
            codes.append(code)
    fp.close()
    return codes

def main():
    fp = codecs.open("classifications.txt", "w", "UTF8")
    codes = read_iso_codes_from_file("iso_codes.txt")
    errors = []
    for index, iso in enumerate(codes):
        try:
            classi = get_classification_from_iso(iso)
            fp.write("%s: %s\n" % (iso, classi))
            print "Got classification for %s: %s (%d/%d)" % (iso, classi, index+1, len(codes) )
            time.sleep(10)
        except Exception as e:
            print "[ERROR] Something went wrong with %s" % iso
            print e
            errors.append(iso)
    fp.close()
    print "Had %d problems in total" % len(errors)
    print "Problematic codes were: %s" % ", ".join(errors)

if __name__ == "__main__":
    main()

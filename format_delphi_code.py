# -*- coding: utf-8 -*-

import copy

def get_text(url=""):
    f = open(url, 'r', encoding="utf-8")
    text = f.read()
    f.close()
    return text


def format_code(text=""):
    tabs_count = 0
    rows = text.split("\n")
    new_text = []
    for i in rows:
        if "end" in i:
            if tabs_count > 0:
                tabs_count -= 1
        pretab = "\t" * tabs_count
        new_text.append(pretab + i.strip())
        if "begin" in i:
            tabs_count += 1
    tabs_count = 0
    rows = copy.deepcopy(new_text)
    new_text = []
    for i in rows:
        if "until" in i:
            if tabs_count > 0:
                tabs_count -= 1
        pretab = "\t" * tabs_count
        new_text.append(pretab + i)
        if "repeat" in i:
            tabs_count += 1
    return "\n".join(new_text)

def save_file(url="", text=""):
    f = open(url, 'w', encoding="utf-8")
    f.write(text)
    f.close()


text = format_code(text=get_text(url="temp_code"))
save_file(url="formated_code.txt", text=text)

'''
if tdict['product']['manufacturer'] is None:
    manufacturer_name = "None"
else:
    if 'name' in tdict['product']['manufacturer']:
        manufacturer_name = tdict['product']['manufacturer']['name']
    else:
        manufacturer_name = "None"

'''



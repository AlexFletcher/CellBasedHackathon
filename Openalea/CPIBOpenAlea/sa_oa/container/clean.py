import os
from os import path

def clean_directory(dir_name):
    list_rep = []
    for filename in os.listdir(dir_name):
        filepath = path.join(dir_name,filename)
        if path.isdir(filepath) :
            list_rep.append(filepath)
        else:
            base, ext = path.splitext(filename)
            if ext == '.pyc':
                os.remove(filepath)
                print filepath

    for rep in list_rep: 
        clean_directory(rep)

if __name__=='__main__' :
    clean_directory(".")

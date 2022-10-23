import os 
def emptyfolder(folder):
    for files in os.listdir(folder):
        path = folder+"/"+files
        if (os.path.isfile(path)==True):
            os.remove(path)

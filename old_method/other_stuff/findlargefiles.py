from pathlib import Path
import os
from utils import listdirs

def find_in_folder(folder):
    names_sizes = []
    for f in folder.iterdir():
        if f.is_dir():
            continue

        name = str(f.absolute())
        size = os.path.getsize(name)
        
        names_sizes.append((name, size))

    return names_sizes
        

def find_in_subfolders(folder):
    names_sizes = []
    for sf in folder.iterdir():
        if not sf.is_dir():
            continue

        names_sizes += find_in_folder(sf)

    return names_sizes


def file_sizes(folder):
    files_here = find_in_folder(folder)
    files_there = find_in_subfolders(folder)

    return files_here + files_there


def run_all(folders):
    r = []
    N = len(folders)
    i = 0
    for f in folders:
        i += 1
        print(f"Running folder {i}/{N}", end="\r")
        r += file_sizes(f)
    return r



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to search")
    parser.add_argument("-z", "--dryrun", action="store_true", help="Do dry-run.")

    args = parser.parse_args()

    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    folders = [f for f in folders if f.is_dir()]
    
    names_sizes = run_all(folders)

    names_sizes = sorted(names_sizes, key=lambda t:t[1])

    for n, s in names_sizes:
        if s / 1e6 < 15:
            continue
        sf = str(round(s / 1e6, 2)) + "MB"
        print(sf, n)
    print("Printed all files above 15 MB")

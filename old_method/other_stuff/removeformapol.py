import os
import sys
from pathlib import Path

def remove_formal(folder, dryrun=False):
    f1 = f"{folder}/formalpol.gpw"
    f2 = f"{folder}/asr.formalpolarization.FAILED"

    if os.path.exists(f1) and os.path.exists(f2):
        if not dryrun:
            print(f"Removing {f1}")
            os.remove(f1)
        else:
            print(f"Would remove {f1}")



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Bilayer folder to cleanup.")
    parser.add_argument("-z", "--dryrun", action="store_true", help="Do dry-run.")
    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]
        

    for folder in folders:
        remove_formal(folder, args.dryrun)

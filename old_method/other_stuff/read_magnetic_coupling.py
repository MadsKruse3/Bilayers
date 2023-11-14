from argparse import ArgumentParser
import numpy as np
import os


parser = ArgumentParser()
add = parser.add_argument
add("-f", "--folder", default=".")

args = parser.parse_args()

numbers = []

for fname in os.listdir(args.folder):
    if "magnetic_bilayer" not in fname: continue
    if ".err" in fname: continue

    number = fname.split(".")[-2]

    numbers.append(int(number))


numbers = sorted(numbers)

fname = f"{args.folder}/magnetic_bilayer.py.{numbers[-1]}.out"

with open(fname, "r") as f:
    txt = f.read()


value = None
for line in txt.split("\n"):
    if line.startswith("eDIFF"):
        value = float(line.split(": ")[-1])
        break

if args.folder == ".":
    mat = os.getcwd().split("/")[-1]
else:
    mat = "/".join(args.folder.split("/")[-2:])

vs = str(value) if value < 0 else "+" + str(value)
print(f"{mat}".ljust(30) + " has eFM - eAFM".ljust(20) + f"= {vs}")

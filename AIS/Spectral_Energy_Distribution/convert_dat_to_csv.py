# # -*- coding: utf-8 -*-
# Nov 3rd 2022
# @author: denis
# This script fix the dat files in the spectral library

import pandas as pd
import numpy as np
import os

for file in os.listdir():
    if file.endswith(".dat"):
        print(file)
        with open(file, "r") as f:
            lines = f.readlines()
        with open("csv//" + file.split(".dat")[0] + ".csv", "w") as f:
            f.write("wavelength (nm),flux (F_lambda)\n")
            for line in lines:
                if line == "":
                    continue
                wv = float(line[:7].replace(" ", "")) / 10
                f.write(str(wv) + "," + line[7:17].replace(" ", "") + "\n")

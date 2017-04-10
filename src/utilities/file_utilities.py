import csv
import os.path

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/utilities")
from error_utilities import error

## Checks to see if a file exists
# @param filename  name of the file to check for existence
def checkFileExists(filename):
  if not os.path.isfile(filename):
    error("The file '" + filename + "' does not exist.")

## Writes data to a CSV file
# @param data  dictionary of lists to output to file
# @param filename  name of the output file
# @param precision  precision of entries in file
def writeCSVFile(data, filename, precision):
  with open(filename, 'wb') as csvfile:
    writer = csv.writer(csvfile)

    # write header row
    key_list = [key for key in data]
    writer.writerow(key_list)

    # create format string
    format_string = "%." + str(precision) + "e"

    # write data
    for i,item in enumerate(data[key_list[0]]):
      row_data = [format_string % (data[key][i]) for key in data]
      writer.writerow(row_data)

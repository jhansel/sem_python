import csv
import numpy as np
import os.path

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

## Reads data from a CSV file
# @param inputfile  name of the file to read
# @return dictionary of lists of each variable
def readCSVData(inputfile):
  first_line_processed = False
  data = dict()
  variable_name_to_index = dict()

  # read data from input file into lists
  with open(inputfile) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
      if (row):
        if (first_line_processed):
          for key in data:
            data[key].append(float(row[variable_name_to_index[key]]))
        else:
          for i,entry in enumerate(row):
            variable_name_to_index[entry] = i
            data[entry] = list()
          first_line_processed = True

  # convert lists to numpy arrays
  for key in data:
    data[key] = np.array(data[key])

  return data

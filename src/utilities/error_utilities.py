import sys
from termcolor import colored

def error(message):
  raise Exception(colored("\nERROR: " + message, "red"))

def inputError(line_number, message):
  error("Input line " + str(line_number) + ": " + message)

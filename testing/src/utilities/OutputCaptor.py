import sys
from cStringIO import StringIO

class OutputCaptor(object):
  def __init__(self):
    # create a backup reference of the standard output
    self.backup = sys.stdout

    # create a new standard output
    sys.stdout = StringIO()

  def getCapturedOutput(self):
    # release the standard output
    out = sys.stdout.getvalue()

    # close the stream and restore the original
    sys.stdout.close()
    sys.stdout = self.backup

    # return the captured output
    return out

if __name__ == "__main__":
  captor = OutputCaptor()
  sys.stdout.write("This line should be third and uppercase\n")
  print "This line should be fourth and uppercase"
  sys.stdout.write("This line should be fifth and uppercase\n")

  out = captor.getCapturedOutput()
  sys.stdout.write("This line should be first and lowercase\n")
  print "This line should be second and lowercase"

  print out.upper()

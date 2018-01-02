#!/usr/bin/env python
import sys
import os
import glob

from optparse import OptionParser


class FofnToSmrtPipeInput:

    """ Construct a smrtpipe input.xml file using a fofn (file of file names) file."""

    def __init__(self, argv):
        self.__parseOptions(argv)

    def __parseOptions(self, argv):
        usage = 'Usage: %prog [--help] [options] bas.fofn > input.xml'
        parser = OptionParser(usage=usage,
                              description=FofnToSmrtPipeInput.__doc__)

        parser.add_option('--jobid',
                          help="Specify a job ID")
        parser.add_option('--jobname',
                          help="Specify a job name")
        parser.add_option('--jobcomment',
                          help="Specify a job comment")

        parser.set_defaults(jobid=None, jobname='',
                            jobcomment='')

        self.options, args = parser.parse_args(argv)

        if len(args) != 2:
            parser.error('Expected 1 argument')
        self.inputFofn = args[1]

    def run(self):
        print '<?xml version="1.0"?>'
        print '<pacbioAnalysisInputs>'

        if self.options.jobid:
            print '<job id="%s">' % self.options.jobid
            print '  <description>'
            print '    <name>%s</name>' % self.options.jobname
            print '    <comment>%s</comment>' % self.options.jobcomment
            print '  </description>'
            print '</job>'

        print '  <dataReferences>'

        infile = open(self.inputFofn, 'r')
        for i, line in enumerate(infile):
            if len(line.strip()) > 0:
                print '    <url ref="run:0000000-%04d">'\
                    '<location>%s</location></url>' % (i, line.strip())
        infile.close()

        print '  </dataReferences>'
        print '</pacbioAnalysisInputs>'

if __name__ == '__main__':
    app = FofnToSmrtPipeInput(sys.argv)
    sys.exit(app.run())

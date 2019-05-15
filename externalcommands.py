import logging
import os
import subprocess


def touch(path):
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    with open(path, 'a'):
        os.utime(path, None)


class ExternalCommand(object):

    def __init__(self, cmd, output_file=None, force_run=False, touch_output=False):

        self.name = os.path.basename(cmd.split(' ')[0]).upper()
        self.force_run = force_run
        self.output_file = output_file
        self.cmd = cmd
        self.touch_output = touch_output

    def run(self):

        if self.output_file is not None and os.path.exists(self.output_file):
            logging.info('SKIPPING COMMAND "{}" - output exists: {}'.format(self.name, self.output_file))
            return 0
        else:
            if self.touch_output:
                touch(self.output_file)

        logging.info('EXEC COMMAND "{}" :\n{}'.format(self.name, self.cmd))
        c = subprocess.call(self.cmd, shell=True)
        logging.info('COMMAND "{}" return code: {}'.format(self.name, c))
        return c

class ExternalCommandOutput(object):

    def __init__(self, cmd):
        #logging.debug(cmd)
        self.cmd = cmd.split(' ')

    def run(self):
        self.proc = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.out, self.err = self.proc.communicate()
        return self.out.decode("utf-8")
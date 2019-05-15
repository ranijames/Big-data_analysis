import logging
import subprocess

__author__ = 'ARJ'


CMD_with_body = 'echo "{}" | mutt -n -s "JOB: {}" {}'
CMD = 'mutt -n -s "JOB: {}" {} < /dev/null'



def send_email(destination, subject, body=None):
    if body is not None:
        #logging.debug('sending email notification with body')
        cmd = CMD_with_body.format(body, subject, destination)
    else:
        #logging.debug('sending email notification')
        cmd = CMD.format(subject, destination)
    c = subprocess.Popen(cmd, shell=True)
    streamdata = c.communicate()[0]
    logging.debug(streamdata)
    logging.debug('Email return code {}'.format(c.returncode))
    return c.returncode


if __name__ == "__main__":  # detects if called from system
    send_email('alvarani@gmail.com', 'Automatic notification test', 'The python script was run from command line')
    send_email('alvarani@gmail.com', 'Automatic notification test without body')

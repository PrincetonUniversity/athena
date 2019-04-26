"""Provides LogPipe class to pipe output from subprocess to a log.
Adapted from https://codereview.stackexchange.com/questions/6567"""
import logging
import threading
import os


class LogPipe(threading.Thread):
    def __init__(self, logger, level):
        """Setup the object with a logger and a loglevel and start the thread"""
        super(LogPipe, self).__init__()
        # threading.Thread.__init__(self)
        self.logger = logging.getLogger(logger)
        self.daemon = False
        self.level = level
        self.fdRead, self.fdWrite = os.pipe()
        self.pipeReader = os.fdopen(self.fdRead)
        self.start()

    def fileno(self):
        """Return the write file descriptor of the pipe"""
        return self.fdWrite

    def run(self):
        """Run the thread, logging everything."""
        for line in iter(self.pipeReader.readline, ''):
            self.logger.log(self.level, line.strip('\n'))
        self.pipeReader.close()

    def close(self):
        """Close the write end of the pipe."""
        os.close(self.fdWrite)

import gzip
import sys
import os
import subprocess


class BasicUtils:

    def read_file(file):
        if (file.endswith('.gz')):
            try:
                fh = gzip.open(file, 'rt')
            except IOError as err:
                print("I/O error: {0}".format(err))
                sys.exit()
        else:
            try:
                fh = open(file, 'r')
            except IOError as err:
                print("I/O error: {0}".format(err))
                sys.exit()
        return (fh)

    def write_file(file):
        if (file.endswith('.gz')):
            try:
                fh = gzip.open(file, 'wt')
            except IOError as err:
                print("I/O error: {0}".format(err))
                sys.exit()
        else:
            try:
                fh = open(file, 'w')
            except IOError as err:
                print("I/O error: {0}".format(err))
                sys.exit()
        return (fh)
    
    def exists(d, keychain):
        # for i in range(0,len(keychain)):
        for k in keychain:
            if (k in d):
                if type(d[k]).__name__ == 'dict':
                     d = d[k]
            else:
                return(False)
        return(True)

    def add(d, keychain , val):
       for k in range(0,len(keychain)):
            if k == len(keychain)-1:
                d[keychain[k]] = val
                return
            if keychain[k] in d:
                if type(d[keychain[k]]).__name__ == 'dict':
                     d = d[keychain[k]]
                else:
                     del d[keychain[k]]
                     d[keychain[k]] = {}
                     d = d[keychain[k]]
            else:
                d[keychain[k]] = {}
                d = d[keychain[k]]

    def make_path(path):
        if not os.path.exists(path):
            os.makedirs(path)

    def get_time():
        now = time.strftime("%c")
        return now

    def run_command(command):
        process = subprocess.Popen(command, shell=True)
        process.wait()

    def runCommand(command, stdout=None, stderr=None):
    	# process = subprocess.Popen(command, shell=True)
    	if stdout is not None:
            ofh = open(stdout, 'w')
    	else:
            ofh = None
    	if stderr is not None:
            efh = open(stderr, 'w')
    	else:
            efh = None

    	process = subprocess.Popen(command, shell=True, stdout=ofh, stderr=efh)
    	process.wait()

"""
This will create and contain the code for maintainng the logs
"""
import os
import logging
import logging.handlers

#create logger
logger = logging.getLogger('CancerSim')
logger.setLevel(logging.DEBUG)

directory =  os.path.dirname(os.path.abspath(__file__))
log_dir = os.path.relpath(os.path.join(directory,'..','logs'))

LOG_FILENAME = 'CancerSim.log'
INFO_FILENAME = 'CancerSim.info'
MAX_BYTES = 2000000

"""
#create file handler
fh = logging.handlers.RotatingFileHandler(
		os.path.join(log_dir,LOG_FILENAME),
		maxBytes=MAX_BYTES,backupCount=5)
fh.setLevel(logging.DEBUG)

#create the info file handler
fh2 = logging.handlers.RotatingFileHandler(
		os.path.join(log_dir,INFO_FILENAME),
		maxBytes=MAX_BYTES,backupCount=5)
fh2.setLevel(logging.INFO)
"""

#create the console logger
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

#create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#fh.setFormatter(formatter)
#fh2.setFormatter(formatter)
ch.setFormatter(formatter)

#logger.addHandler(fh)
#logger.addHandler(fh2)
logger.addHandler(ch)

logger.info('Created the logging instance...')

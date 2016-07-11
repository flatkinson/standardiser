import logging

####################################################################################################

# Defaults...

level   = 'INFO'
fmt     = '[%(asctime)s %(name)s %(levelname)s] %(message)s'
datefmt = '%d/%m/%y %H:%M:%S'

####################################################################################################

def run(name, level=level, fmt=fmt, datefmt=datefmt):

	logger = logging.getLogger(name)

	handler = logging.StreamHandler()

	handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))

	logger.addHandler(handler)

	logger.setLevel(level)

	logger.propagate = False

	return logger

# run

####################################################################################################
# End
####################################################################################################

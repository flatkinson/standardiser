# Turn on verbose logging in ipython: %run -i debug.py
import logging
reload(logging)
logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")

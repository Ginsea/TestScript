import logging

def setlog(logfile=None):
    if logfile:
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)

        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        logger.addHandler(fh)
        logger.addHandler(ch)
    else:
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        ch.setFormatter(formatter)

        logger.addHandler(ch)
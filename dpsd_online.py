#!/usr/bin/env python

import sys, time, os, logging
import dpsd_run, dicxml
from aug_sfutils import journal, getlastshot, SFREAD

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSDon')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

firstshot = False
# On a non-shotday the loop is terminated by the crontab final time (7 p.m.)
while not firstshot:
    firstshot = journal.anyshot()
    if not firstshot:
        time.sleep(60)

dpsd_dir = os.path.dirname(os.path.realpath(__file__))
xml_d = dicxml.xml2dict('%s/xml/shot.xml' %dpsd_dir)['main']
setup_d = dicxml.xml2val_dic(xml_d)


if journal.anyshot():

    hour = 0
    while (hour < 19):

        loctime = time.localtime(time.time())
        hour = loctime[3]
        logger.info('Waiting for next shot')
        logger.info('ctrl+c to terminate the script')

        try:
            nshot = getlastshot()
            logger.info('Last shot: %d' %nshot)
        except:
#            raise ValueError('Problems reading last shot number')
            logger.warning('Problems reading last shot number')
            time.sleep(60)
            continue
 
        try:
            nsp = SFREAD(nshot, 'NSP')
            if nsp.status:
                logger.warning('NSP shotfile for shot %d exists already, writing no new shotfile' %nshot)
            else:
                setup_d['io']['Shots'] = nshot 
                dpsd_run.DPSD(setup_d)
                setup_d['io']['Shots'] = nshot - 1
                dpsd_run.DPSD(setup_d)
        except:
            logger.error('Problems writing shotfile for #%d' %nshot)

        time.sleep(60)
 
else:

    logger.info('Today is no shotday')

sys.exit()


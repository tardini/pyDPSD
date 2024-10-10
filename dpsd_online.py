#!/usr/bin/env python

import sys, time, os, logging, json
import dpsd_run
import aug_sfutils as sf

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSDon')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

dpsd_dir = os.path.dirname(os.path.realpath(__file__))
f_json = '%s/settings/shot.json' %dpsd_dir
with open(f_json) as fjson:
    setup_d = json.load(fjson)

hour = 0
while (hour < 19):
    loctime = time.localtime(time.time())
    hour = loctime[3]
    logger.info('Waiting for next shot')
    logger.info('ctrl+c to terminate the script')

    lastshot = sf.wait4shot()
    try:
        nsp = sf.SFREAD(lastshot, 'NSP')
        if nsp.status:
            logger.warning('NSP shotfile for shot %d exists already, writing no new shotfile' %nshot)
        else:
            setup_d['io']['Shots'] = lastshot 
            dpsd_run.DPSD(setup_d)
            setup_d['io']['Shots'] = lastshot - 1
            dpsd_run.DPSD(setup_d)
    except:
        logger.error('Problems writing shotfile for #%d' %nshot)

    time.sleep(60)


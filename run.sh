#!/bin/bash
echo "create event data"
./event.exe ES.21833.FULL.DAQ_TEST.20190824163728.001.dat ES.21833.FULL.DAQ_TEST.20190824163728.001.dat.event.root

echo "create status data"
./status.exe ES.21833.FULL.DAQ_TEST.20190824163728.001.dat ES.21833.FULL.DAQ_TEST.20190824163728.001.dat.status.root

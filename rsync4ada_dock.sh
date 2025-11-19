#!/bin/bash

# Create local directories if they don't exist
mkdir -p /data/work/dock/virtual_screening/results/
mkdir -p /data/work/dock/virtual_screening/docking_output/
mkdir -p /data/work/dock/virtual_screening/input/receptors/
mkdir -p /data/work/dock/analisis/
mkdir -p /data/work/dock/vs_crystal/

rsync -avvz --progress --update --times -e "ssh"  \
    --exclude='*.npz'\
    --exclude='*.xml'\
    --exclude='*.dlg'\
    jaeohshin@syntax.kias.re.kr:/store/jaeohshin/work/dock/ \
    /data/work/dock/
exit

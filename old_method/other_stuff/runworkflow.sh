#!bin/bash

source ~/.bashrc
source ~/standardvenv/bin/activate

if [ `mq ls ~/stacked/tree/ -sqr | wc -l` -lt 500 ]
then
    python3 -W ignore -m myqueue workflow ~/bilayerworkflow/revisedworkflow.py ~/stacked/tree/AB2/VI2/VI2-d480b343c982
fi

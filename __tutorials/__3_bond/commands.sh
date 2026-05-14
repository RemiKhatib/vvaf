#!/bin/bash

rm *dat *xyz

awk -f pos.awk > pos.xyz
awk -f vel.awk pos.xyz > vel.xyz

../../vvaf input_vvaf


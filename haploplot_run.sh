#!/bin/bash

read -p "Please state sex of input sample (m/f) " RESP
if [ "$RESP" = "m" ]; then
  sh ./scripts/xy_haploplot.sh
else
  sh ./scripts/xx_haploplot.sh
fi
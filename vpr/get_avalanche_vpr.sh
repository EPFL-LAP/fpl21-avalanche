#!/bin/bash

wget https://github.com/verilog-to-routing/vtr-verilog-to-routing/archive/refs/tags/v8.0.0.tar.gz
tar -zxvf v8.0.0.tar.gz
patch -s -p1 < vtr8_avalanche.patch

#!/bin/env bash

../defcell.py -s exx -a 10.0 -o out.pdb test.pse
diff out.pdb test_exx10.pdb

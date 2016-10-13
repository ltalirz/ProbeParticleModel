#!/bin/bash

echo "test for the PP-STM code:"
python PPSTM/gen_LJFF.py -i crazy_mol.xyz
python PPSTM/rel_scan.py --pos
python PPSTM/plot_AFM_res.py --pos --df --save_df
python PPSTM/dIdV_example.py
echo "Now all things made, removing not neccessary folders and *.npy files"
rm -r ./Q0.00K0.50
rm *.npy
echo "Done, bye bye"



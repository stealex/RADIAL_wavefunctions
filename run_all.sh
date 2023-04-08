#!/bin/bash

#nuclei=("78Kr" "96Ru" "102Pd" "106Cd" "124Xe" "130Ba" "136Ce")
nuclei=("78Kr" "96Ru" "106Cd" "124Xe" "130Ba" "136Ce" "50Cr" "58Ni" "62Zn" "74Se" "84Sr" "92Mo" "102Pd" "112Sn" "120Te" "144Sm" "156Dy" "162Er" "168Yb" "174Hf" "184Os" "190Pt" "36Ar" "40Ca" "54Fe" "108Cd" "126Xe" "132Ba" "138Ce" "152Gd" "158Dy" "164Er" "180W" "196Hg")
for nuc in ${nuclei[@]}; do
  rm -rf "${nuc}_2EC"
  python dhfs.py "dhfs_${nuc}_2EC.ini"
done
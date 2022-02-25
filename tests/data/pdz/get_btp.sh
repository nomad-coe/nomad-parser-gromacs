

sys="2awx"

/sw/linux/gromacs/5.1.2/bin/gmxdump -s $sys.tpr &> $sys.dumped.txt

/data/bee7/rudzinski/soft_backup/BOCS_Nov_2019/BOCS/force-matching/build/bin/translator -f $sys.dumped.txt -o $sys.btp



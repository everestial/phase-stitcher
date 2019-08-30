python3 phase-Stitcher.py --nt 1 --input example_01/haplotype_file01.txt --output outdir --mat MA605 --pat Sp21 --f1Sample ms02g --culLH maxSum --lods 3

# Namespace(chr='', culLH='maxSum', f1Sample='ms02g', hapStats='no', 
# input='example_01/haplotype_file01.txt', lods='3', mat='MA605', nt='1',
#  outPatMatID='pat,mat', output='outdir', pat='Sp21')


python3 phase-Stitcher.py --nt 2 --input example_01/haplotype_file01.txt --output tests/outdir/egout2 --mat MA605,MA622 --pat Sp21,Sp164 --f1Sample ms02g --outPatMatID Sp,My --culLH maxSum --lods 3

# Namespace(chr='', culLH='maxSum', f1Sample='ms02g', hapStats='no', 
# input='example_01/haplotype_file01.txt', lods='3', mat='MA605,MA622', nt='2',
# outPatMatID='Sp,My', output='tests/outdir/egout2', pat='Sp21,Sp164')

python3 phase-Stitcher.py --nt 1 --input example_02/haplotype_file02.txt --output tests/outdir/egout3 --mat pre:MA,Nc --pat pre:Sp --f1Sample ms02g --culLH maxPd --lods 25
# Namespace(chr='', culLH='maxPd', f1Sample='ms02g', hapStats='no',
# input='example_02/haplotype_file02.txt', lods='25', mat='pre:MA,Nc', nt='1', 
# outPatMatID='pat,mat', output='tests/outdir/egout3', pat='pre:Sp')

python3 phase-Stitcher.py --nt 1 --input example_02/haplotype_file02.txt --output tests/outdir/egout4 --mat pre:MA,Nc --pat pre:Sp --f1Sample ms02g --culLH maxPd --lods 25 --hapStats yes --outPatMatID Sp,My

# Namespace(chr='', culLH='maxPd', f1Sample='ms02g', hapStats='yes', 
# input='example_02/haplotype_file02.txt', lods='25', mat='pre:MA,Nc', nt='1', 
# outPatMatID='Sp,My', output='tests/outdir/egout4', pat='pre:Sp')
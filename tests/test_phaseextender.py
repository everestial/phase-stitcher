import filecmp
import os
import shutil
import tempfile
import pathlib
import difflib
import webbrowser

import pytest
from collections import namedtuple


from phaser import phase_stich

soi_meta_cls = namedtuple('soi_meta_cls','soif1, soimom, soidad, dad_id, mom_id')

# input_file = 'tests/inputs/haplotype_file01.txt'
# soif1 = 'ms02g'
# soimom = [('MA605:PI', 'MA605:PG_al')]
# soidad = [('Sp21:PI', 'Sp21:PG_al')]
# dad_id = 'pat_hap'
# mom_id =  'mat_hap'
# outputdir = 'tests/temp1'
# chr_list = ''
# nt = 1
# lods_cut_off = 3
# max_is = 'max sum'
# maxed_as = '+'
# hapstats  = 'no'


def test_first_example():
    with tempfile.TemporaryDirectory() as temp_dir1:
        soi_meta = soi_meta_cls('ms02g', [('MA605:PI', 'MA605:PG_al')], [('Sp21:PI', 'Sp21:PG_al')], 'pat_hap', 'mat_hap')
        args1 = ('tests/inputs/haplotype_file01.txt', soi_meta, temp_dir1, '', 1, 3, '+', 'no')
        phase_stich(*args1)
        assert is_same_dir('tests/outdir/egout1/', temp_dir1) 


def test_second_example():
    with tempfile.TemporaryDirectory() as temp_dir1:
        soi_meta = soi_meta_cls('ms02g', [('MA605:PI', 'MA605:PG_al'), ('MA622:PI', 'MA622:PG_al')], [('Sp21:PI', 'Sp21:PG_al'), ('Sp164:PI', 'Sp164:PG_al')], 'Sp_hap', 'My_hap',)
        args2 = ('tests/inputs/haplotype_file01.txt', soi_meta, temp_dir1, '', 2, 3, '+', 'no')
        phase_stich(*args2)
        assert is_same_dir('tests/outdir/egout2', temp_dir1 ) 


def test_third_example():
    with tempfile.TemporaryDirectory() as temp_dir1:
        soi_meta = soi_meta_cls('ms02g', [('MA605:PI', 'MA605:PG_al'), ('MA625:PI', 'MA625:PG_al'), ('MA611:PI', 'MA611:PG_al'), ('MA629:PI', 'MA629:PG_al'), ('MA622:PI', 'MA622:PG_al'), ('Ncm8:PI', 'Ncm8:PG_al')], [('Sp76:PI', 'Sp76:PG_al'), ('Sp154:PI', 'Sp154:PG_al'), ('Sp164:PI', 'Sp164:PG_al'), ('Sp3:PI', 'Sp3:PG_al'), ('SpNor33:PI', 'SpNor33:PG_al'), ('Sp21:PI', 'Sp21:PG_al')], 'pat_hap', 'mat_hap')
        args3 = ('tests/inputs/haplotype_file02.txt', soi_meta, temp_dir1, '', 1, 25, '*', 'no')
        phase_stich(*args3)
        assert is_same_dir('tests/outdir/egout3', temp_dir1 ) 


def test_fourth_example():
    with tempfile.TemporaryDirectory() as temp_dir1:
        soi_meta = soi_meta_cls('ms02g', [('MA611:PI', 'MA611:PG_al'), ('MA625:PI', 'MA625:PG_al'), ('Ncm8:PI', 'Ncm8:PG_al'), ('MA622:PI', 'MA622:PG_al'), ('MA629:PI', 'MA629:PG_al'), ('MA605:PI', 'MA605:PG_al')], [('Sp154:PI', 'Sp154:PG_al'), ('Sp164:PI', 'Sp164:PG_al'), ('Sp76:PI', 'Sp76:PG_al'), ('Sp21:PI', 'Sp21:PG_al'), ('SpNor33:PI', 'SpNor33:PG_al'), ('Sp3:PI', 'Sp3:PG_al')], 'Sp_hap', 'My_hap')
        args4 = ('tests/inputs/haplotype_file02.txt', soi_meta, temp_dir1, '', 1, 25, '*', 'yes')
        phase_stich(*args4)
        assert is_same_dir('tests/outdir/egout4', temp_dir1 ) 


def is_same_dir(dir1, dir2):
    """
    Compare two directory trees content.
    Return False if they differ, True is they are the same.
    """
    compared = filecmp.dircmp(dir1, dir2)
    if (compared.left_only or compared.right_only or compared.diff_files
        or compared.funny_files):
        print(compared.left_only)
        print(compared.right_only)
        print(compared.diff_files)

        # NOTE Uncomment below lines if you want to see diff between two files
        # it is only when test cases fail
        # for f in compared.diff_files:
        #     diff_html(dir1+'/'+f, dir2+'/'+f)
        return False
    for subdir in compared.common_dirs:
        if not is_same_dir(os.path.join(dir1, subdir), os.path.join(dir2, subdir)):
            return False
    return True 


def diff_html(filepath1, filepath2):
    """
    This function is useful while debugging tests failures.
    It compares two files line by line in browser 
    """

    fromlines = pathlib.Path(filepath1).read_text().splitlines()
    tolines = pathlib.Path(filepath2).read_text().splitlines()

    diff = difflib.HtmlDiff(tabsize=4,).make_file(fromlines, tolines)
    # is_same_dir(file1, file2)
    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
        url = 'file://' + f.name
        f.write(diff)
    webbrowser.open_new_tab(url)
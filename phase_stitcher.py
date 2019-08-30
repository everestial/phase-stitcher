import sys
import argparse

from args_builder import get_args
from args_extractor import args_to_values
from utils import print_authorship
from phaser import phase_stich

def main():
    parser = argparse.ArgumentParser()
    args = get_args(parser)
    parsed_args = args_to_values(args)
    print_authorship()
    phase_stich(*parsed_args)

if __name__ == "__main__":
    main()
    


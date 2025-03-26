#!/usr/bin/env python3
"""
Homological Threading - Main Interface

This script provides a unified interface to access analysis and plotting functionality.
Usage:
    python main.py analysis [options]
    python main.py plot [options]
"""

import argparse
import sys
import pathlib
from importlib import import_module


def run_analysis(args):
    """Run analysis script with given arguments"""
    # Import the analysis module
    sys.path.append(str(pathlib.Path(__file__).resolve().parent))
    analysis_module = import_module("scripts.analysis")
    
    # Replace the first argument with the script name
    sys.argv = [str(pathlib.Path("scripts/analysis.py"))] + args
    
    # Call the main function of the analysis module
    analysis_module.main()


def run_plot(args):
    """Run plot script with given arguments"""
    # Import the plot module
    sys.path.append(str(pathlib.Path(__file__).resolve().parent))
    plot_module = import_module("scripts.plot")
    
    # Replace the first argument with the script name
    sys.argv = [str(pathlib.Path("scripts/plot.py"))] + args
    
    # Call the main function of the plot module
    plot_module.main(plot_module.get_args())


def main():
    # Create the parser
    parser = argparse.ArgumentParser(
        description='Homological Threading Interface',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py analysis pd -i data/*.data -o output/
  python main.py analysis betti -i output/*.h5 -o output/ -f results.h5
  python main.py analysis num_threading -i output/*.h5 -o output/ -f results.h5
  python main.py plot pd -i output/*.h5
  python main.py plot betti -i output/results.h5
  python main.py plot num_threading -i output/results.h5 --type pdf
"""
    )
    
    # Add subparsers for analysis and plot
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    subparsers.required = True
    
    # Analysis command
    analysis_parser = subparsers.add_parser('analysis', help='Run analysis functions')
    analysis_parser.add_argument('args', nargs=argparse.REMAINDER, help='Arguments for analysis script')
    
    # Plot command
    plot_parser = subparsers.add_parser('plot', help='Run plotting functions')
    plot_parser.add_argument('args', nargs=argparse.REMAINDER, help='Arguments for plot script')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Run the appropriate command
    if args.command == 'analysis':
        run_analysis(args.args)
    elif args.command == 'plot':
        run_plot(args.args)


if __name__ == '__main__':
    main() 
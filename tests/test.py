import matplotlib.pyplot as plt
import sys
import pathlib
import time
import argparse
import numpy as np
import os
import h5py

sys.path.append(str(pathlib.Path(__file__).resolve().parent.parent / "src"))

import homological_threading as ht


def test(filename, output):
    """
    Run homological threading test on a LAMMPS data file.
    
    args:
    filename: str
        Input LAMMPS data file.
    output: str
        Output HDF5 file.
    """
    print(f"Testing homological threading on {filename}")
    
    # Initialize HomologicalThreading
    time_start = time.time()
    pds = ht.HomologicalThreading()
    coords = pds.read_lmpdata(filename)
    print(f"Read {len(coords)} chains from {filename}")
    
    # Compute persistence diagrams for individual chains
    time_pd_i_start = time.time()
    pds.pd_i.compute(coords, dim=1, mp=False)
    time_pd_i_end = time.time()
    print(f"Elapsed time for computing pd_i: {time_pd_i_end - time_pd_i_start:.2f} seconds")
    
    # Compute persistence diagrams for pairs of chains
    time_pd_i_cup_j_start = time.time()
    pds.pd_i_cup_j.compute(coords, dim=1, mp=True)
    time_pd_i_cup_j_end = time.time()
    print(f"Elapsed time for computing pd_i_cup_j: {time_pd_i_cup_j_end - time_pd_i_cup_j_start:.2f} seconds")
    
    # Compute threading
    time_threading_start = time.time()
    pds.threading.compute(pds.pd_i.pd, pds.pd_i_cup_j.pd)
    time_threading_end = time.time()
    print(f"Elapsed time for computing threading: {time_threading_end - time_threading_start:.2f} seconds")
    
    # Get number of chains
    nchains = pds.pd_i.pd.shape[0]
    print(f"Number of chains: {nchains}")
    
    # Process pd_i
    all_pds_i = pds.pd_i.pd.reshape(-1, 2)
    print(f"pd_i shape: {pds.pd_i.pd.shape}, reshaped: {all_pds_i.shape}")
    
    # Process pd_i_cup_j
    all_pds_i_cup_j = pds.pd_i_cup_j.pd.reshape(-1, 2)
    print(f"pd_i_cup_j shape: {pds.pd_i_cup_j.pd.shape}, reshaped: {all_pds_i_cup_j.shape}")
    
    # Process threading
    all_pds_threading = pds.threading.pd.reshape(-1, 2)
    print(f"threading shape: {pds.threading.pd.shape}, reshaped: {all_pds_threading.shape}")
    
    # Filter out invalid points (birth time <= 0 or death time <= birth time)
    all_pds_i = all_pds_i[all_pds_i[:, 0] > 0]
    all_pds_i = all_pds_i[all_pds_i[:, 1] > all_pds_i[:, 0]]
    
    all_pds_i_cup_j = all_pds_i_cup_j[all_pds_i_cup_j[:, 0] > 0]
    all_pds_i_cup_j = all_pds_i_cup_j[all_pds_i_cup_j[:, 1] > all_pds_i_cup_j[:, 0]]
    
    all_pds_threading = all_pds_threading[all_pds_threading[:, 0] > 0]
    all_pds_threading = all_pds_threading[all_pds_threading[:, 1] > all_pds_threading[:, 0]]
    
    print(f"Valid points - pd_i: {len(all_pds_i)}, pd_i_cup_j: {len(all_pds_i_cup_j)}, threading: {len(all_pds_threading)}")
    
    # Save results to HDF5 file using the class method
    pds.to_hdf5(output)
    print(f"Results saved to {output}")
    
    # Calculate Betti numbers using the class methods
    alphas_i, betti_i = pds.pd_i.betti()
    alphas_i_cup_j, betti_i_cup_j = pds.pd_i_cup_j.betti()
    alphas_threading, betti_threading = pds.threading.betti()
    
    return {
        "pd_i": (alphas_i, betti_i),
        "pd_i_cup_j": (alphas_i_cup_j, betti_i_cup_j),
        "threading": (alphas_threading, betti_threading)
    }


def validate_hdf5(file_path):
    """
    Validate the HDF5 file structure.
    
    args:
    file_path: str
        Path to the HDF5 file.
    
    returns:
    bool: True if valid, False otherwise.
    """
    try:
        # Create a HomologicalThreading instance and load from HDF5
        pds = ht.HomologicalThreading()
        pds.from_hdf5(file_path)
        
        # Check if all required data is loaded
        if pds.pd_i.pd is None:
            print("Missing dataset: pd_i")
            return False
            
        if pds.pd_i_cup_j.pd is None:
            print("Missing dataset: pd_i_cup_j")
            return False
            
        if pds.threading.pd is None:
            print("Missing dataset: threading")
            return False
        
        # Check shapes
        nchains = pds.pd_i.pd.shape[0]
        
        if pds.pd_i_cup_j.pd.shape[:2] != (nchains, nchains):
            print(f"Shape mismatch in pd_i_cup_j: expected first dimensions ({nchains}, {nchains}), got {pds.pd_i_cup_j.pd.shape[:2]}")
            return False
        
        if pds.threading.pd.shape[:2] != (nchains, nchains):
            print(f"Shape mismatch in threading: expected first dimensions ({nchains}, {nchains}), got {pds.threading.pd.shape[:2]}")
            return False
        
        print(f"HDF5 file validation successful:")
        print(f"  Number of chains: {nchains}")
        print(f"  pd_i shape: {pds.pd_i.pd.shape}")
        print(f"  pd_i_cup_j shape: {pds.pd_i_cup_j.pd.shape}")
        print(f"  threading shape: {pds.threading.pd.shape}")
        
        # Print metadata if available
        if pds.metadata:
            print("  Metadata:")
            for key, value in pds.metadata.items():
                if value is not None:
                    print(f"    {key}: {value}")
        
        return True
    except Exception as e:
        print(f"Error validating HDF5 file: {e}")
        return False

def batch_test(input_dir, output_dir, pattern="*.data"):
    """
    Run tests on multiple input files.
    
    args:
    input_dir: str
        Directory containing input files.
    output_dir: str
        Directory for output files.
    pattern: str
        Glob pattern for input files.
    """
    import glob
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find input files
    input_files = glob.glob(os.path.join(input_dir, pattern))
    print(f"Found {len(input_files)} input files matching pattern '{pattern}'")
    
    results = {}
    for input_file in input_files:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(output_dir, f"{base_name}.h5")
        
        print(f"\n{'='*80}")
        print(f"Processing {input_file}")
        print(f"{'='*80}")
        
        try:
            result = test(
                input_file, 
                output_file, 
                plot=False, 
                save_data=True
            )
            
            # Validate output file
            valid = validate_hdf5(output_file)
            results[base_name] = {
                "success": True,
                "valid_hdf5": valid,
                "data": result
            }
        except Exception as e:
            print(f"Error processing {input_file}: {e}")
            results[base_name] = {
                "success": False,
                "error": str(e)
            }
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"Summary")
    print(f"{'='*80}")
    
    success_count = sum(1 for r in results.values() if r["success"])
    print(f"Processed {len(results)} files, {success_count} successful, {len(results) - success_count} failed")
    
    for name, result in results.items():
        status = "SUCCESS" if result["success"] else f"FAILED: {result['error']}"
        print(f"  {name}: {status}")
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test homological threading")
    parser.add_argument("--input", "-i", type=str, help="Input LAMMPS data file or directory")
    parser.add_argument("--output", "-o", type=str, help="Output HDF5 file or directory")
    parser.add_argument("--batch", "-b", action="store_true", help="Batch process all files in input directory")
    parser.add_argument("--pattern", "-p", type=str, default="*.data", help="File pattern for batch processing")
    
    args = parser.parse_args()
    
    proj_dir = pathlib.Path(__file__).resolve().parent.parent
    
    # Default values if not specified
    if args.input is None:
        args.input = str(proj_dir / "data/N10M100.data")
    if args.output is None:
        args.output = str(proj_dir / "data/ht.h5")
    
    # Batch processing
    if args.batch:
        if not os.path.isdir(args.input):
            print(f"Error: Input must be a directory for batch processing")
            sys.exit(1)
        
        if not os.path.isdir(args.output):
            os.makedirs(args.output, exist_ok=True)
        
        batch_test(args.input, args.output, args.pattern)
    else:
        # Single file processing
        test(
            args.input, 
            args.output, 
        )
        
        # Validate output file
        if os.path.exists(args.output):
            validate_hdf5(args.output)

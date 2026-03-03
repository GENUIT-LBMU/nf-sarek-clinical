#!/usr/bin/env python3
"""
Script to run HRProfiler with command line arguments
Usage: python HRD.py --snv-dir <path> --cnv-dir <path> --output-dir <path> [options]
"""

import argparse
import sys
from HRProfiler.scripts import HRProfiler as HR

def main():
    parser = argparse.ArgumentParser(
        description='Run HRProfiler for Homologous Recombination Deficiency analysis'
    )
    
    parser.add_argument(
        '--snv-dir',
        required=True,
        help='Directory with filtered VCF files'
    )
    
    parser.add_argument(
        '--cnv-dir',
        required=True,
        help='Directory with segment files'
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to save results'
    )
    
    parser.add_argument(
        '--cnv-file-type',
        default='ASCAT',
        choices=['ASCAT', 'SEQUENZA', 'FACETS', 'PURPLE'],
        help='CNV file type (default: ASCAT)'
    )
    
    parser.add_argument(
        '--organ',
        default='BREAST',
        choices=['BREAST', 'OVARIAN'],
        help='Organ type for prediction (default: BREAST)'
    )
    
    parser.add_argument(
        '--genome',
        default='GRCh38',
        help='Reference genome version (default: GRCh38)'
    )
    
    parser.add_argument(
        '--hrd-threshold',
        type=float,
        default=0.5,
        help='HRD probability threshold (default: 0.5)'
    )
    
    parser.add_argument(
        '--nreplicates',
        type=int,
        default=20,
        help='Number of bootstrap replicates (default: 20)'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("HRProfiler - Homologous Recombination Deficiency Analysis")
    print("=" * 60)
    print(f"\nParameters:")
    print(f"  SNV Directory:    {args.snv_dir}")
    print(f"  CNV Directory:    {args.cnv_dir}")
    print(f"  CNV File Type:    {args.cnv_file_type}")
    print(f"  Output Directory: {args.output_dir}")
    print(f"  Organ Type:       {args.organ}")
    print(f"  Genome:           {args.genome}")
    print(f"  HRD Threshold:    {args.hrd_threshold}")
    print(f"  N Replicates:     {args.nreplicates}")
    print("\nStarting analysis...\n")
    
    try:
        HR.HRProfiler(
            genome=args.genome,
            exome=True,
            INDELS_DIR=None,
            SNV_DIR=args.snv_dir,
            CNV_DIR=args.cnv_dir,
            RESULT_DIR=args.output_dir,
            cnv_file_type=args.cnv_file_type,
            bootstrap=False,
            nreplicates=args.nreplicates,
            normalize=True,
            hrd_prob_thresh=args.hrd_threshold,
            plot_predictions=True,
            organ=args.organ
        )
        
        print("\n" + "=" * 60)
        print("✓ HRProfiler completed successfully!")
        print("=" * 60)
        print(f"\nResults saved to: {args.output_dir}\n")
        
    except Exception as e:
        print(f"\n✗ Error running HRProfiler: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

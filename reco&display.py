#!/usr/bin/env python3

import argparse
import subprocess

def run_o2_mch_reco_workflow():
    # Construct the command for o2-mch-reco-workflow
    reco_workflow_command = [
        "o2-mch-reco-workflow",
        "-b",
        "--condition-remap", "http://localhost:6464=MCH/Calib/RejectList",
        #"--configKeyValues", "MCHStatusMap.useRejectList=true;MCHDigitFilter.statusMask=2",
        "|",
        "o2-qc",
        "--config", "json://./qc-mch-clusters.json",
        "--local-batch=QC.root",
        "|",
        "o2-dpl-run",
        "-b"
    ]

    try:
        # Run o2-mch-reco-workflow
        subprocess.run(" ".join(reco_workflow_command), shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}")

def run_o2_mch_map_mch(green, norm_per_area, rootfile_left, rootfile_right):
    # Construct the command for o2-mch-map_mch
    map_mch_command = [
        "o2-mch-map_mch",
        "--green" if green else "",
        "--normperarea" if norm_per_area else "",
        "--rootfileleft", rootfile_left,
        "--rootfileright", rootfile_right
    ]

    # Remove empty strings from the command (resulting from optional arguments)
    map_mch_command = [arg for arg in map_mch_command if arg]

    try:
        # Run o2-mch-map_mch
        subprocess.run(map_mch_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}")

def main():
    parser = argparse.ArgumentParser(description='Run o2-mch-reco-workflow or o2-mch-map_mch')
    parser.add_argument('command', choices=['display', 'reconstruction'], default='display', nargs='?',
                        help='Choose a command to run (default: display)')

    args = parser.parse_args()

    if args.command == 'reconstruction':
        run_o2_mch_reco_workflow()
    elif args.command == 'display':
        run_o2_mch_map_mch(green=True, norm_per_area=True, rootfile_left="DATA_QC.root", rootfile_right="100mil.root")

if __name__ == "__main__":
    main()

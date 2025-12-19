#!/usr/bin/env python3

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter tandem repeats (TRs) based on the tools used to detect them.")
    parser.add_argument("-t", "--trf", required=True, help="Path to the TRF input file.")
    parser.add_argument("-r", "--rm", required=True, help="Path to the RepeatMasker output file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TRF filtered file.")
    return parser.parse_args()


def parse_gff_line(line):
    """Parse a GFF line into a dictionary."""
    fields = line.strip().split("\t")
    return {
        "seqid": fields[0],
        "type": fields[2],
        "start": int(fields[3]),
        "end": int(fields[4]),
    }


def check_overlap(tr_entry, rm_entry):
    """Check if two regions overlap"""
    start1 = min(tr_entry["start"], rm_entry["start"])
    start2 = max(tr_entry["start"], rm_entry["start"])
    end1 = min(tr_entry["end"], rm_entry["end"])
    end2 = max(tr_entry["end"], rm_entry["end"])

    return start1 <= end2 and end1 >= start2


def main(): 
    args = parse_arguments()

    # Read the RepeatMasker file
    with open(args.rm, 'r') as file:
        rm_dict = {}
        for line in file:
            if not line.startswith("#"):  # Ignore header lines
                rm_entry = parse_gff_line(line)
                if rm_entry["type"] == "tandem_repeat":  # Check if it's a tandem repeat
                    rm_dict.setdefault(rm_entry["seqid"], []).append(rm_entry)

    # Filter TRF lines by matching with RepeatMasker
    with open(args.trf, 'r') as trf_file, open(args.output, 'w') as output_file:
        removed_count = 0
        for line in trf_file:
            if line.startswith("#"):
                # Write comment lines directly to the output
                output_file.write(line)

            else:
                # Process TRF line
                trf_entry = parse_gff_line(line)
                seqid = trf_entry["seqid"]
                overlapping = False

                if seqid in rm_dict:
                    for rm_entry in rm_dict[seqid]:
                        if check_overlap(trf_entry, rm_entry):
                            removed_count += 1
                            overlapping = True
                            break
                
                if not overlapping:
                    output_file.write(line)

    print("Filtering complete. Results saved to", args.output, "with", removed_count, "TRs removed.")

if __name__ == "__main__":
    main()

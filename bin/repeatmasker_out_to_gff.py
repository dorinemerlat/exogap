#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert RepeatMasker .out to GFF with 9 columns"
    )
    p.add_argument('-r', '--repeatmasker', required=True, help='RepeatMasker .out file')
    p.add_argument('-m', '--mchelper', required=False, help='mchelper output (optional mapping file)')
    p.add_argument('-o', '--output', required=True, help='Output GFF file')
    return p.parse_args()


def load_mchelper_map(path: str):
    """Load MCHelper TSV as a dict key -> row-dict.
    Index by multiple identifiers so lookups succeed:
      - id (e.g. RHYIM_TE196)
      - mchelper_id or MCHelper_id (e.g. RHYIM_TE196_inc)
      - rm2_id when present
    Also expects annotation columns: rm_code, rm2_class, rm2_order, rm2_superfamily, rm2_family, rm2_subfamily, rm2_subfamily2.
    Missing values are handled later as NA.
    """
    if not path:
        return {}
    mp: dict[str, dict] = {}
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
            header = None
            for raw in fh:
                line = raw.rstrip('\n')
                if not line:
                    continue
                parts = line.split('\t')
                # first non-empty line that has many columns -> header
                if header is None:
                    header = parts
                    continue
                # pad to header length
                if len(parts) < len(header):
                    parts = parts + [''] * (len(header) - len(parts))
                row = {header[i]: parts[i] for i in range(len(header))}
                # add using multiple possible keys
                keys = set()
                for k in ('id', 'mchelper_id', 'MCHelper_id', 'rm2_id'):
                    v = row.get(k)
                    if v:
                        keys.add(v)
                for k in keys:
                    mp[k] = row
    except Exception:
        # best-effort: ignore parsing errors
        pass
    return mp


def classify_type(classfamily: str) -> str:
    """Map RepeatMasker class/family to GFF type field."""
    if classfamily is None:
        return 'Transposable_elements'
    cf = classfamily.strip()
    if cf.lower() == 'unknown':
        return 'Transposable_elements'
    if cf.startswith('Simple_repeat') or cf == 'Simple_repeat':
        return 'Simple_repeat'
    if cf.startswith('Low_complexity') or cf == 'Low_complexity':
        return 'Low_complexity'
    return 'Transposable_elements'


def parse_rm_line(tokens):
    """Parse tokens from a RepeatMasker .out line (split by whitespace).
    Returns dict with needed fields or None if malformed.
    Expected columns (simplified):
      0: score
      1: perc div.
      2: perc del.
      3: perc ins.
      4: sequence (seqid)
      5: begin
      6: end
      7: (left)
      8: strand (+ or C)
      9: repeat name
      10: class/family
    """
    try:
        score = tokens[0]
        seqid = tokens[4]
        start = int(tokens[5])
        end = int(tokens[6])
        strand_tok = tokens[8]
        repeat_name = tokens[9]
        classfamily = tokens[10] if len(tokens) > 10 else None
        # positions in repeat (begin, end) may include parentheses for C strand; strip them
        rep_begin = tokens[11] if len(tokens) > 11 else ''
        rep_end = tokens[12] if len(tokens) > 12 else ''
        try:
            rep_begin_i = int(rep_begin.strip('()')) if rep_begin else None
        except Exception:
            rep_begin_i = None
        try:
            rep_end_i = int(rep_end.strip('()')) if rep_end else None
        except Exception:
            rep_end_i = None

        strand = '+' if strand_tok == '+' else ('-' if strand_tok.upper() == 'C' else '.')
        return {
            'score': score,
            'seqid': seqid,
            'start': start,
            'end': end,
            'strand': strand,
            'repeat_name': repeat_name,
            'classfamily': classfamily,
            'rep_begin': rep_begin_i,
            'rep_end': rep_end_i,
        }
    except Exception:
        return None

def main():
    args = parse_args()
    in_path = Path(args.repeatmasker)
    if not in_path.exists():
        print(f"Error: input .out not found: {in_path}", file=sys.stderr)
        sys.exit(2)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mcp = load_mchelper_map(args.mchelper)

    written = 0
    idx = 0
    with in_path.open('r', encoding='utf-8', errors='ignore') as fh, out_path.open('w', encoding='utf-8') as outfh:
        for line in fh:
            # skip header or empty lines
            if not line.strip() or line.startswith('SW') or line.startswith('score'):
                continue
            tokens = line.strip().split()
            # some .out files have leading spaces; ensure at least 11 tokens
            if len(tokens) < 11:
                continue
            parsed = parse_rm_line(tokens)
            if not parsed:
                continue

            idx += 1
            seqid = parsed['seqid']
            start = parsed['start']
            end = parsed['end']
            score = parsed['score']
            strand = parsed['strand']
            repeat_name = parsed['repeat_name']

            # class/family from .out (fallback) — we keep for GFF type
            classfamily = parsed['classfamily']
            gff_type = classify_type(classfamily)

            # phase is '.' for repeats
            phase = '.'

            # Build attribute field according to requested spec
            total_length = end - start + 1
            # repeat type and code
            repeatmasker_code = (classfamily or 'NA')
            repeatmasker_code_norm = repeatmasker_code.capitalize() if repeatmasker_code != 'NA' else 'NA'
            # determine repeat_type
            if classfamily is None or classfamily.lower() == 'unknown':
                repeat_type = 'Unknown'
                repeat_order = 'NA'
                repeat_superfamily = 'NA'
            elif classfamily.startswith('Simple_repeat'):
                repeat_type = 'Simple_repeat'
                repeat_order = 'NA'
                repeat_superfamily = 'NA'
            elif classfamily.startswith('Low_complexity'):
                repeat_type = 'Low_complexity'
                repeat_order = 'NA'
                repeat_superfamily = 'NA'
            else:
                repeat_type = 'Transposable_elements'
                # parse order/superfamily from classfamily like LTR/Gypsy
                parts_cf = classfamily.split('/')
                repeat_order = parts_cf[0] if parts_cf else 'NA'
                repeat_superfamily = parts_cf[1] if len(parts_cf) > 1 else 'NA'

            matching_repeat = repeat_name
            # Build ID as RE_<prefix>_<idx>
            prefix = seqid.split('_')[0] if '_' in seqid else seqid
            repeat_id = f"RE_{prefix}_{idx}"

            # If an MCHelper row exists for this repeat id, use it to fill fields
            mrow = mcp.get(matching_repeat, {}) if isinstance(mcp, dict) else {}
            def val(row, key, default='NA'):
                v = row.get(key, '') if isinstance(row, dict) else ''
                return v if v not in (None, '', 'NA', 'na') else default

            rm2_class = val(mrow, 'rm2_class', 'NA')
            rm2_order = val(mrow, 'rm2_order', 'NA')
            rm2_super = val(mrow, 'rm2_superfamily', 'NA')
            rm2_family = val(mrow, 'rm2_family', 'NA')
            rm2_subfam = val(mrow, 'rm2_subfamily', 'NA')
            rm2_subfam2 = val(mrow, 'rm2_subfamily2', 'NA')
            rm_code = val(mrow, 'rm_code', 'NA')
            # prefer mchelper row values; otherwise keep computed values
            out_repeat_class = rm2_class if rm2_class != 'NA' else 'NA'
            out_repeat_order = rm2_order if rm2_order != 'NA' else repeat_order
            out_repeat_super = rm2_super if rm2_super != 'NA' else repeat_superfamily
            out_repeat_family = rm2_family if rm2_family != 'NA' else 'NA'
            out_repeat_subfam = rm2_subfam if rm2_subfam != 'NA' else 'NA'
            out_repeat_subfam2 = rm2_subfam2 if rm2_subfam2 != 'NA' else 'NA'
            out_rm_code = rm_code if rm_code != 'NA' else repeatmasker_code_norm

            attrs = [
                f"ID={repeat_id}",
                f"total_length={total_length}",
                f"matching_repeat={matching_repeat}",
                f"repeat_type={repeat_type}",
                f"repeat_class={out_repeat_class}",
                f"repeat_order={out_repeat_order}",
                f"repeat_superfamily={out_repeat_super}",
                f"repeat_family={out_repeat_family}",
                f"repeat_subfamily={out_repeat_subfam}",
                f"repeat_subfamily2={out_repeat_subfam2}",
                f"repeatmasker_code={out_rm_code}",
                f"position_in_repeat_begin={parsed['rep_begin'] if parsed['rep_begin'] is not None else 'NA'}",
                f"position_in_repeat_end={parsed['rep_end'] if parsed['rep_end'] is not None else 'NA'}",
            ]
            # ensure trailing semicolon
            attr_str = ';'.join(attrs) + ';'

            # columns: seqid, source, type, start, end, score, strand, phase, attributes
            outfh.write('\t'.join([
                seqid,
                'RepeatMasker',
                gff_type,
                str(start),
                str(end),
                str(score),
                strand,
                phase,
                attr_str,
            ]) + '\n')
            written += 1

    print(f"✅ Wrote {written} records to {out_path}")


if __name__ == '__main__':
    main()

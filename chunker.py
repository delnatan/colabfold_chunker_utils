import argparse
from pathlib import Path


def parse_fasta(filename):
    sequences = {}

    with open(filename, "r") as f:
        current_id = None
        current_sequence = ""

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = current_sequence
                current_id = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_id is not None:
            sequences[current_id] = current_sequence

    return sequences


def save_chunks(sequence, id_str, segment_length, overlap_length, output_dir,
                line_width=80):

    seqlength = len(sequence)
    chunk_step = segment_length - overlap_length

    if segment_length >= seqlength:
        s = 0
        e = seqlength
        segment_sequence = sequence[s:e]
        subsegments = [
            segment_sequence[i:i+line_width]
            for i in range(0, len(segment_sequence), line_width)
        ]
        with open(output_dir / ("%s.fa" % (id_str)), "wt") as f:
            f.write(">%s_%d-%d\n" % (id_str, s+1, e))
            for segment in subsegments:
                f.write(segment)
                f.write("\n")
            f.write("\n")

        # do early termination of function (return nothing)
        return

    nchunks = seqlength // chunk_step + 1

    outfilename = "%s_segments.fa" % (id_str)

    # define end of sequence index (iterate until we reach the end)
    e = 0
    i = 0
    with open(output_dir / outfilename, "wt") as f:
        while e < seqlength:
            s = i * (segment_length - overlap_length)
            e = min(s + segment_length, seqlength)

            segment_sequence = sequence[s:e]

            # break sequence up to conform to 'line_width'
            subsegments = [
                segment_sequence[i:i+line_width]
                for i in range(0, len(segment_sequence), line_width)
            ]

            f.write(">%s_seg%d_%d-%d\n" % (id_str, i+1, s+1, e))
            for subsegment in subsegments:
                f.write(subsegment)
                f.write("\n")

            f.write("\n")
            
            # increment chunk counter
            i += 1


def main():
    desc_str = r'''
    Divides up long protein sequence into overlapping chunks

                      overlap_length
                     __|__
                    |--------------------|
    |---------------------|
        segment_length

    '''
    parser = argparse.ArgumentParser(
        prog="chunker v1",
        description=desc_str,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        argument_default=argparse.SUPPRESS
    )

    parser.add_argument("fasta_file", type=str, help="input fasta file")
    parser.add_argument("-segment_length", default=1400, type=int,
                        help="segment length (default: %(default)d)")
    parser.add_argument("-overlap_length", default=200, type=int,
                        help="overlap length (default: %(default)d)")
    parser.add_argument("--output_dir", default=".", type=str,
                        help="output directory (default: %(default)s)")

    args = parser.parse_args()

    sequences = parse_fasta(args.fasta_file)
    output_dir = Path(args.output_dir)

    for seq_id, sequence in sequences.items():
        save_chunks(sequence, seq_id, args.segment_length,
                    args.overlap_length, output_dir)


if __name__ == "__main__":

    main()

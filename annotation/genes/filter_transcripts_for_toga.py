#!/usr/bin/env python3
"""Filter transcripts for TOGA input.

Remove:
- noncoding trasncripts
- transcripts with NMD targets (optional)
- broken transcripts (adjustable threshold)

Requires:
- Annotation File
- Isoforms files (will be generated, if not provided)
"""
import argparse
import sys
import os
import re
from collections import Counter
from collections import defaultdict
from datetime import datetime as dt
import networkx as nx
from twobitreader import TwoBitFile

__author__ = "Bogdan Kirilenko, 2021"
__version__ = "0.2"

# TODO: fusion transcript detection
# TODO: output formatting, test bed and genePred outputs
# TODO: add another condition to identical CDS: picket transcript MUST NOT be a NMD

BED = "bed"
GENEPRED = "genepred"
ORIGINAL = "original"
BED_SIZE = 12
GENEPRED_SIZE = 15
ISOFORMS_FILE_COLS = 2

ALLOWED_CHARSET = "a-zA-Z0-9._-"
ALLOWED_CHARSET_RE = rf"[^{ALLOWED_CHARSET}]"

SUPPORTED_IN_FMT = {BED, GENEPRED}
SUPPORTED_OUT_FMT = {ORIGINAL, BED, GENEPRED}
PLUS = "+"
MINUS = "-"
DOT = "."
STRANDS = {PLUS, MINUS, DOT}

START = "ATG"
STOPS = {"TAG", "TGA", "TAA"}
TGA = "TGA"
ACCEPTOR_SITES = {"AG"}
DONOR_SITES = {
    "GT",
    "GC",
}
U12_ACCEPTOR_SITES = {"AT"}
U12_DONOR_SITES = {"AC"}
NMD_DISTANCE = 55
MAX_NUM_SKIPPED_GENES_TO_CONSOLE = 50

complement = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "G",
    "n": "n",
}


class Transcript:
    """Internal transcript representation."""

    def __init__(self, line, fmt=BED):
        """Initiate the class reading a line."""
        self.raw_line = line  # keep the raw line
        self.in_fmt = fmt
        trans_raw_data = {}  # define here to supress linter
        if fmt.lower() == BED:
            # depending on the input format, we get a bit different data
            # parsing bed-12 format
            trans_raw_data = self.__parse_bed(line)
        elif fmt.lower() == GENEPRED:
            # reading genePredEx (15 columns)
            trans_raw_data = self.__parse_genepred(line)
        # TODO: reading standard genePred (12 col)
        else:  # not supported: raise an error
            raise ValueError(f"Transcript class: not supported format {fmt}")
        # trans_raw_data dict now has all data obtained from the line

        # these parameters are mandatory:
        # chromosome id, transcription start and end,
        # transcript identifier
        self.chrom = trans_raw_data["chrom"]
        self.chrom_start = trans_raw_data["chrom_start"]
        self.chrom_end = trans_raw_data["chrom_end"]
        self.id = trans_raw_data["name"]
        self.gene_locus_size = abs(self.chrom_end - self.chrom_start)
        # if non-unique and we filter non-uniq transcripts out then
        # write here what replaced this transcript
        self.non_uniq_replaced_with = None

        # not really necessary, bed or genePred-specific fields
        self.bed_score = trans_raw_data.get("bed_score", 0)
        self.item_rgb = trans_raw_data.get("item_rgb", "0,0,0")
        self.name_2 = trans_raw_data.get("name_2", "None")
        self.cds_start_stat = trans_raw_data.get("cds_start_stat", "unknown")
        self.cds_end_stat = trans_raw_data.get("cds_end_stat", "unknown")
        # may be necessary for bed->genePred conversion
        self.exon_frames = trans_raw_data.get("exon_frames", None)

        # necessary but can be omitted (which is not recommended)
        self.strand = trans_raw_data.get("strand", "+")
        # if translation start is not provided: use transcription start
        self.trans_start = trans_raw_data.get("thick_start", self.chrom_start)
        self.trans_end = trans_raw_data.get("thick_end", self.chrom_end)

        # exons data: to be read later
        self.exons_num = trans_raw_data["block_count"]
        self.exon_starts = trans_raw_data["block_starts"]
        self.exon_ends = trans_raw_data["block_ends"]
        self.exon_sizes = trans_raw_data["block_sizes"]

        # get CDS part; set default values (non-coding)
        # overwrite if some cds exons detected
        self.exon_utr_or_cds = []  # U - fully UTR, C - fully CDS exon, UC - both
        self.cds_exon_starts = []
        self.cds_exon_ends = []
        self.cds_exon_sizes = []
        self.cds_exons_num = 0
        self.utr_left_len = 0
        self.utr_rigth_len = 0
        self.__get_cds_exons()  # fills self.cds_*** variables
        self.is_nmd__no_seq = False  # default val
        self.__check_nmd_no_seq()
        self.tot_cds_len = sum(self.cds_exon_sizes)
        # infer some necessary data
        self.frame_is_correct = self.tot_cds_len % 3 == 0
        self.is_coding = self.tot_cds_len > 0 or self.strand == "."
        self.intron_sizes = []
        self.min_intron_size = 0  # if 0 -> single exon gene
        self.__get_intron_sizes()  # fill introns data
        if self.exon_frames is None:
            # exon phases were not read: probably input is a bed file
            self.__get_exon_frames()
        # if None -> was not set
        # False -> not, True -> yes
        self.corresponding_gene = None

        # flags: whether something is wrong
        self.starts_with_atg = None
        self.ends_with_stop = None
        self.is_nmd = None
        # flag: if False, then 3 variables above cannot be used
        self.seq_was_checked = False
        self.u12_introns = []
        self.broken_introns = []
        self.inframe_stop_codons = []
        # for selenocysteine-filtes
        self.all_stops_are_tga = None
        self.__is_intact = None

        # ralation to other isoforms
        self.is_longest_intact_cds_in_gene = None
        self.is_long_enough = None  # is not shorter than X% than the longest
        self.is_principal = None
        # maybe ids of isoforms of the same gene?

        # user explicitly asked to add this transcript:
        self.must_include = False

        # decision about this transcript
        self.selected = None
        # need to highlight this scenario:
        self.selected_to_avoid_skipping_gene = False

    def __parse_bed(self, line):
        """Parse bed-12-formattedline."""
        line_data = line.rstrip().split("\t")
        # bed could be 12, 9, 4 or whatever
        # for now we read only the bed 12
        if len(line_data) != BED_SIZE:
            raise ValueError(
                f"Bed file expected to have {BED_SIZE} fields, got:\n{line}"
            )

        try:  # raises ValueError if the problem is known
            # then the except block catches this and adds the general message
            # that shows the corrupted line itself
            chrom = line_data[0]
            chrom_start = int(line_data[1])
            chrom_end = int(line_data[2])
            if chrom_end < chrom_start:
                raise ValueError("chromEnd < chromStart")
            name = line_data[3]
            bed_score = int(line_data[4])
            strand = line_data[5]
            if strand not in STRANDS:
                raise ValueError(f"Unknown strand {strand}")
            thick_start = int(line_data[6])
            thick_end = int(line_data[7])
            if thick_end < thick_start:
                raise ValueError("thickEnd < thickStart")
            item_rgb = line_data[8]
            block_count = int(line_data[9])
            block_sizes = [int(x) for x in line_data[10].split(",") if x != ""]
            block_starts = [int(x) for x in line_data[11].split(",") if x != ""]
            if block_count != len(block_sizes) or block_count != len(block_starts):
                raise ValueError("blockCount field is wrong")
            block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
        except ValueError as verr:
            sys.stderr.write(str(verr))
            sys.stderr.write("\n")
            sys.stderr.write(f"Corrupted line:\n{line}")
            sys.exit(1)

        ret = {
            "chrom": chrom,
            "chrom_start": chrom_start,
            "chrom_end": chrom_end,
            "name": name,
            "bed_score": bed_score,
            "strand": strand,
            "thick_start": thick_start,
            "thick_end": thick_end,
            "item_rgb": item_rgb,
            "block_count": block_count,
            "block_sizes": block_sizes,
            "block_starts": block_starts,
            "block_ends": block_ends,
        }
        return ret

    def __parse_genepred(self, line):
        """Parse genePred line."""
        line_data = line.rstrip().split("\t")
        if len(line_data) != GENEPRED_SIZE:
            raise ValueError(
                f"genePred file expected to have {GENEPRED_SIZE} fields, got:\n{line}"
            )

        try:
            name = line_data[0]
            chrom = line_data[1]
            strand = line_data[2]
            if strand not in STRANDS:
                raise ValueError(f"Unknown strand {strand}")
            chrom_start = int(line_data[3])
            chrom_end = int(line_data[4])
            if chrom_end < chrom_start:
                raise ValueError("chromEnd < chromStart")
            thick_start = int(line_data[5])
            thick_end = int(line_data[6])
            if thick_end < thick_start:
                raise ValueError("thickEnd < thickStart")
            block_count = int(line_data[7])
            block_starts_ = [int(x) for x in line_data[8].split(",") if x != ""]
            block_ends_ = [int(x) for x in line_data[9].split(",") if x != ""]
            # it was initially written for the bed format, which provides relative exon coords
            # so here we convert them to relatoves
            block_starts = [x - chrom_start for x in block_starts_]
            block_ends = [x - chrom_start for x in block_ends_]
            block_sizes = [abs(x[1] - x[0]) for x in zip(block_starts, block_ends)]
            bed_score = int(line_data[10])
            name_2 = line_data[11]
            cds_start_stat = line_data[12]
            cds_end_stat = line_data[13]
            exon_frames = [int(x) for x in line_data[14].split(",") if x != ""]
        except ValueError as verr:
            sys.stderr.write(str(verr))
            sys.stderr.write("\n")
            sys.stderr.write(f"Corrupted line:\n{line}")
            sys.exit(1)

        ret = {
            "chrom": chrom,
            "chrom_start": chrom_start,
            "chrom_end": chrom_end,
            "name": name,
            "bed_score": bed_score,
            "strand": strand,
            "thick_start": thick_start,
            "thick_end": thick_end,
            "block_count": block_count,
            "block_sizes": block_sizes,
            "block_starts": block_starts,
            "block_ends": block_ends,
            "name_2": name_2,
            "cds_start_stat": cds_start_stat,
            "cds_end_stat": cds_end_stat,
            "exon_frames": exon_frames,
        }
        return ret

    @staticmethod
    def __count_U_in_a_row(lst):
        """Get a number of UTR exons in a row from the beginning of the list."""
        ret = 0
        for x in lst:
            if x != "U":
                break
            ret += 1
        return ret

    def __get_exon_frames(self):
        """Compute exon frames acording to UCSC genePredExt format."""
        if len(self.cds_exon_sizes) == 0:
            # no CDS: all exons are non-coding
            self.exon_frames = [-1 for _ in range(self.exons_num)]
            return
        utr_exons_left = self.__count_U_in_a_row(self.exon_utr_or_cds)
        utr_exons_rigth = self.__count_U_in_a_row(self.exon_utr_or_cds[::-1])
        if len(self.cds_exon_sizes) == 1:
            # nothing to compute for a single CDS-exon transcript
            self.exon_frames = (
                [-1 for _ in range(utr_exons_left)]
                + [0]
                + [-1 for _ in range(utr_exons_rigth)]
            )
        phases = [-1 for _ in range(utr_exons_left)]
        left_pointer = 0
        right_pointer = 0
        cds_phases = []
        # if +: read exon sizes from the beginning to the end
        # if -: from end to the beginning
        iter_over_ = (
            self.cds_exon_sizes if self.strand == PLUS else self.cds_exon_sizes[::-1]
        )
        for size in iter_over_:
            cds_phases.append(right_pointer)
            right_pointer = (size - left_pointer) % 3
            left_pointer = 3 - right_pointer
        # need to invert phases if - phase
        cds_phases = cds_phases if self.strand == PLUS else cds_phases[::-1]
        phases.extend(cds_phases)
        phases.extend([-1 for _ in range(utr_exons_rigth)])
        self.exon_frames = [x for x in phases]

    def __get_intron_sizes(self):
        """Set intron sizes."""
        # We ignore UTR introns!
        if self.cds_exons_num < 2:
            # if so: nothing we can do
            return
        # make a list like
        # exon 1 start exon 1 end exon 2 start exon 2 end etc... remove 1st and last elements
        # split into pairs, the result will be:
        # (exon 1 end, exon 2 start), (exon 2 end, exon 3 start), etc..
        # previously we required at least 2 exons
        # then split list into pairs -> each pair will be (intron_start, intron_end)
        breaking_points = flatten(zip(self.cds_exon_starts, self.cds_exon_ends))[1:-1]
        intron_coords = parts(breaking_points, 2)
        self.intron_sizes = [abs(x[1] - x[0]) for x in intron_coords]
        self.min_intron_size = min(self.intron_sizes)

    def __check_nmd_no_seq(self):
        """Check whether this is an NMD.

        Method uses only the annotation data,
        no seq data.
        """
        # get relative start of the last exon (including UTRs)
        if self.exons_num < 2:  # there is no way
            # for now we assume that 0 exons is also possible
            return
        # get relative coordinate of translation end and the last junction
        if self.strand == PLUS:
            transl_end = self.utr_left_len + sum(self.cds_exon_sizes)
            last_junction_coord = sum(self.exon_sizes[:-1])
            distance = last_junction_coord - transl_end
        else:
            transl_end = self.utr_left_len
            last_junction_coord = self.exon_sizes[0]
            distance = transl_end - last_junction_coord
        # if distance is negative -> stop in the last exon
        if distance > NMD_DISTANCE:
            self.is_nmd__no_seq = True

    def __get_cds_exons(self):
        """Trim UTRs, get CDS exons."""
        if self.trans_start == self.trans_end:
            # nothing to do, left default values
            return
        block_abs_starts = [
            self.exon_starts[i] + self.chrom_start for i in range(self.exons_num)
        ]
        block_abs_ends = [
            self.exon_ends[i] + self.chrom_start for i in range(self.exons_num)
        ]

        for block_num in range(self.exons_num):
            # go block-by-block
            block_start = block_abs_starts[block_num]
            block_end = block_abs_ends[block_num]
            block_size = abs(block_start - block_end)

            # skip the block if it is entirely UTR
            if block_end <= self.trans_start:
                self.utr_left_len += block_size
                self.exon_utr_or_cds.append("U")
                continue
            elif block_start >= self.trans_end:
                self.utr_rigth_len += block_size
                self.exon_utr_or_cds.append("U")
                continue

            # block is either CDS entirely or
            # check block start
            if block_start >= self.trans_start:
                block_new_start = block_start
                self.exon_utr_or_cds.append("C")
            else:
                # start intersects thickStart, part of exon is UTR
                block_new_start = self.trans_start
                utr_piece_ = abs(self.trans_start - block_start)
                self.exon_utr_or_cds.append("UC")
                self.utr_left_len += utr_piece_
            # check block end
            if block_end <= self.trans_end:
                block_new_end = block_end
                self.exon_utr_or_cds.append("C")
            else:
                # partly UTR exon, end intersects thickEnd
                block_new_end = self.trans_end
                utr_piece_ = abs(block_end - self.trans_end)
                self.exon_utr_or_cds.append("UC")
                self.utr_rigth_len += utr_piece_
            # save blocks with updated coordinates
            # also convert them back to relative coordinates with - thick_start
            # after the update thick_start/End are equal to chrom_start/End
            self.cds_exon_starts.append(block_new_start - self.trans_start)
            self.cds_exon_ends.append(block_new_end - self.trans_start)
        self.cds_exons_num = len(self.cds_exon_starts)
        self.cds_exon_sizes = [
            self.cds_exon_ends[i] - self.cds_exon_starts[i]
            for i in range(self.cds_exons_num)
        ]

    def get_abs_exons_coords(self):
        """Return a list of exons - absolute coordinates."""
        ret = []
        for i in range(self.cds_exons_num):
            # get absolute coordinates for each exon
            exon_abs_start = self.cds_exon_starts[i] + self.trans_start
            exon_abs_end = exon_abs_start + self.cds_exon_sizes[i]
            exon_tup = (self.chrom, self.strand, exon_abs_start, exon_abs_end)
            ret.append(exon_tup)
        return ret

    def seq_analysis(self, translated_part):
        """Analyze codon sequences."""
        exon_seqs = []
        self.seq_was_checked = True
        for ex_coord in zip(self.cds_exon_starts, self.cds_exon_ends):
            exon_seq = translated_part[ex_coord[0] : ex_coord[1]]
            exon_seqs.append(exon_seq)
        if self.strand == MINUS:
            exon_seqs = [invert_complement(x) for x in exon_seqs][::-1]
        breaking_points = flatten(zip(self.cds_exon_starts, self.cds_exon_ends))[1:-1]
        intron_coords = parts(breaking_points, 2)
        splice_sites = []
        for intron in intron_coords:
            intron_start = translated_part[intron[0] : intron[0] + 2].upper()
            intron_end = translated_part[intron[1] - 2 : intron[1]].upper()
            splice_sites.append((intron_start, intron_end))
        # if negative strand: revert splice sites
        if self.strand == MINUS:
            splice_sites = [
                (invert_complement(x[1]), invert_complement(x[0])) for x in splice_sites
            ][::-1]
        for num, sst in enumerate(splice_sites):
            don_ = sst[0]
            acc_ = sst[1]
            intron_id = f"{num + 1}_{don_}->{acc_}"
            don_canon = don_ in DONOR_SITES
            don_u12 = don_ in U12_DONOR_SITES
            acc_canon = acc_ in ACCEPTOR_SITES
            acc_u12 = acc_ in U12_ACCEPTOR_SITES

            if don_canon is True and acc_canon is True:
                pass  # absolutely normal canonical splice site
            elif don_u12 is True and acc_u12 is True:
                # normal U12 intron
                self.u12_introns.append(intron_id)
            else:
                # something different
                self.broken_introns.append(intron_id)

        coding_sequence = "".join(exon_seqs)
        codons = parts(coding_sequence, 3)
        if len(codons) == 0:
            self.non_coding = True  # must not happen actually
            return

        self.starts_with_atg = codons[0] == START
        self.ends_with_stop = codons[-1] in STOPS
        # get numbers of all stop codons
        num_of_tga = 0
        for num, codon in enumerate(codons[:-1], 2):
            if codon in STOPS:
                self.inframe_stop_codons.append(num)
            if codon == TGA:
                # if we allow selenocysteines:
                # need to track TGA codons separately
                num_of_tga += 1
        if len(self.inframe_stop_codons) > 0:
            if len(self.inframe_stop_codons) == num_of_tga:
                self.all_stops_are_tga = True
            else:
                self.all_stops_are_tga = False
        # otherwise this field is not applicable

    def is_intact(self, allow_tga=False):
        """Return True if intact, False otherwise."""
        if self.seq_was_checked:
            # we can use data inferred from sequence
            inframe_stop = len(self.inframe_stop_codons) > 0
            if self.all_stops_are_tga is True and allow_tga is True:
                inframe_stop = False

            broken_frame = (
                self.starts_with_atg is False
                or self.ends_with_stop is False
                or inframe_stop is True
            )
        else:
            broken_frame = False
        not_intact = (
            broken_frame is True
            or self.is_nmd__no_seq is True
            or self.frame_is_correct is False
            or self.tot_cds_len == 0
        )
        self.__is_intact = not not_intact
        return self.__is_intact

    def output_bed(self, trim_utr=False):
        """Output as a bed line.

        According to:
        http://genome.ucsc.edu/FAQ/FAQformat#format1
        """
        # bed-specific fields:
        # bed_score: default 0
        # item_rgb: default 0,0,0
        exon_sizes_field = ",".join([str(x) for x in self.exon_sizes])
        exon_starts_field = ",".join([str(x) for x in self.exon_starts])
        fields = (
            self.chrom,
            self.chrom_start,
            self.chrom_end,
            self.id,
            self.bed_score,
            self.strand,
            self.trans_start,
            self.trans_end,
            self.item_rgb,
            self.exons_num,
            exon_sizes_field,
            exon_starts_field,
        )
        bed_string = "\t".join([str(x) for x in fields])
        return bed_string

    def output_genepred(self):
        """Output as a gene pred extended line.

        According to:
        http://genome.ucsc.edu/FAQ/FAQformat#format9
        """
        # need to correct coordinates of exons
        # from relative to absolute
        exon_starts_field = (
            ",".join([str(x + self.chrom_start) for x in self.exon_starts]) + ","
        )
        exon_ends_field = (
            ",".join([str(x + self.chrom_start) for x in self.exon_ends]) + ","
        )
        frames_field = ",".join([str(x) for x in self.exon_frames]) + ","
        # fill the necessary fields
        # genePredEx specific:
        # cds_start_stat
        # cds_end_stat
        # name_2
        # exon_frames
        fields = (
            self.id,
            self.chrom,
            self.strand,
            self.chrom_start,
            self.chrom_end,
            self.trans_start,
            self.trans_end,
            self.exons_num,
            exon_starts_field,
            exon_ends_field,
            self.bed_score,
            self.name_2,
            self.cds_start_stat,
            self.cds_end_stat,
            frames_field,
        )
        gene_pred_string = "\t".join([str(x) for x in fields])
        return gene_pred_string

    def __repr__(self):
        """Describe what do we output when printing an instance."""
        out = ["<Transcript object"]
        for k, v in self.__dict__.items():
            out.append(f"{k}: {v}")
        out.append(">")
        return "\n".join(out)

    def __hash__(self):
        """Use transcript id to hash.

        We assume that they are unique.
        """
        return hash(self.id)

    def get_cds_hash(self):
        """Get unique CDS hash."""
        _starts = tuple(self.cds_exon_starts)
        _sizes = tuple(self.cds_exon_sizes)
        to_hash = (_starts, _sizes)
        return hash(to_hash)


def print_stderr(msg, end="\n"):
    """Print to stderr stream."""
    sys.stderr.write(str(msg))
    sys.stderr.write(end)


def flatten(lst):
    """Flat list out of list of lists."""
    return [item for sublist in lst for item in sublist]


def parts(lst, n=3):
    """Split an iterable into parts with size n."""
    return [lst[i : i + n] for i in iter(range(0, len(lst), n))]


def read_transcripts(input_file, input_fmt):
    """Read the input file."""
    out = {}
    f = open(input_file, "r")
    print(f"Reading input files, expected format: {input_fmt}")
    count_unique_transcripts = Counter()
    allowed_re = re.compile(ALLOWED_CHARSET_RE).search
    transcripts_broken_ids = []
    non_coding = 0
    lines_num = 0

    for line in f:
        if line == "\n":
            continue
        trans = Transcript(line, input_fmt)
        count_unique_transcripts[trans.id] += 1
        corr_name = not bool(allowed_re(trans.id))
        if corr_name is False:
            transcripts_broken_ids.append(trans.id)
        lines_num += 1
        if trans.is_coding is False:
            # do not include transcripts without annotated CDS
            # under no circumstances
            non_coding += 1
            continue
        out[trans.id] = trans
    f.close()
    print(f"Processed {lines_num} lines")
    print(f"Skipped {non_coding} non-coding transcripts (CDS length = 0)")
    print(f"Loaded {len(out.keys())} transcripts to filter")

    if lines_num == 0:  # Quit if the input file is empty, write a message about
        print_stderr(f"Warning! {input_file} is empty!")
        sys.exit(0)
    # check whether there are non-uniq transcripts
    # for now non-unique transcripts are not allowed
    # anyway, toga crashes if finds any of them
    if len(transcripts_broken_ids) > 0:
        print_stderr(
            "Error! Some of the transcripts have not allowed characters in the ID:"
        )
        for t in transcripts_broken_ids:
            print_stderr(t)
        print_stderr(f"Only the following characters are allowed: {ALLOWED_CHARSET}")
        sys.exit(1)

    non_uniq = {k for k, v in count_unique_transcripts.items() if v > 1}
    if len(non_uniq) == 0:
        return out
    # there are non-unique transcripts, stop
    print_stderr("Detected non unique transcripts!")
    for elem in non_uniq:
        print_stderr(elem)
    print_stderr("Abort")
    sys.exit(1)


def get_exons_data(trans_dct):
    """Extract exons data to infer isoforms."""
    exons_list = []
    exon_id_to_trans = {}
    exon_counter = 0  # each exon gets an unique ID

    for trans, trans_obj in trans_dct.items():
        t_exons = trans_obj.get_abs_exons_coords()
        for ex in t_exons:
            exon_counter += 1
            _chrom = ex[0]
            _strand = ex[1]
            _start = ex[2]
            _end = ex[3]
            enumed_ex = (exon_counter, _chrom, _strand, _start, _end)
            exons_list.append(enumed_ex)
            exon_id_to_trans[exon_counter] = trans

    return exons_list, exon_id_to_trans


def split_exons_in_chr_dir(exons_list):
    """Make a (chr, strand): [exons] dict."""
    chr_dir_exons_not_sorted = defaultdict(list)
    for exon in exons_list:
        # just rearrange a data a bit
        # go exon-by-exon
        exon_id = exon[0]
        chrom = exon[1]
        strand = exon[2]
        start = exon[3]
        end = exon[4]
        # split into two elements: key (chr & direction) and value (exon range + exon ID)
        chr_dir = (chrom, strand)
        exon_reduced = (exon_id, start, end)
        # add this to default dict
        chr_dir_exons_not_sorted[chr_dir].append(exon_reduced)
    # sort exons in each chr_dir track -> for efficient algorithm
    # use start (field 1) as key
    chr_dir_exons = {
        k: sorted(v, key=lambda x: x[1]) for k, v in chr_dir_exons_not_sorted.items()
    }
    return chr_dir_exons


def intersect_ranges(region_1, region_2):
    """Return TRUE if ranges intersect."""
    inter_size = min(region_1[1], region_2[1]) - max(region_1[0], region_2[0])
    return True if inter_size > 0 else False


def intersect_exons(chr_dir_exons, exon_id_to_transcript):
    """Create graph of connected exons.

    We will use a graph where nodes are transcript IDs.
    Nodes are connected if they have a pair of intersected exons.
    """
    G = nx.Graph()  # init the graph
    for exons in chr_dir_exons.values():
        # this is the same chrom and direction now
        # add nodes now: to avoid missing transcripts
        # that don't intersect anything
        exon_ids = [e[0] for e in exons]
        transcripts = set(exon_id_to_transcript[x] for x in exon_ids)
        G.add_nodes_from(transcripts)

        if len(exons) == 1:
            # only one exon: nothing to intersect
            # gene id is on graph now: we can just continue
            continue
        exons_num = len(exons)
        for i in range(exons_num - 1):
            # nearly N^2 algorithm (a nested loop)
            # but actually faster
            # pick exon i
            exon_i = exons[i]
            exon_i_id = exon_i[0]
            exon_i_start = exon_i[1]
            exon_i_end = exon_i[2]
            # get range and corresponding transcript
            i_range = (exon_i_start, exon_i_end)
            i_trans = exon_id_to_transcript[exon_i_id]
            # then let's look at next exons
            for j in range(i + 1, exons_num):
                # exon j exists for sure (first range is up to exons_num - 1)
                # num of exons is > 1
                exon_j = exons[j]
                exon_j_id = exon_j[0]
                exon_j_start = exon_j[1]
                exon_j_end = exon_j[2]
                if exon_j_start > exon_i_end:
                    # if exon_j start is > exon_i end then we can break
                    # they are sorted: for all next exons start would be > exon_i end
                    # and if so, they do not intersect for sure
                    break
                # let's check that they intersect
                j_range = (exon_j_start, exon_j_end)
                j_trans = exon_id_to_transcript[exon_j_id]  # j exon transcript
                they_intersect = intersect_ranges(i_range, j_range)
                if they_intersect is False:
                    # don't intersect: nothing to do
                    continue
                # if we are in this branch: exons intersect
                # we can add an edge connecting their corresponding transcripts
                G.add_edge(i_trans, j_trans)
    return G


def get_graph_components(graph):
    """Split graph in connected components.

    Different networkx versions have different methods
    to split a graph into connected components.
    Linter will definitely report an error here, depending on the nx version.
    """
    version = nx.__version__.split(".")
    fst = int(version[0])
    snd = int(version[1]) if len(version) > 1 else 9
    num_ = float(f"{fst}.{snd}")
    if num_ < 2.4:
        # old nx version: most likely linter complains here
        graph_components = list(nx.connected_component_subgraphs(graph))
    else:
        graph_components = [graph.subgraph(c) for c in nx.connected_components(graph)]
    return graph_components


def parse_components(conn_components):
    """Parse graph connected components."""
    gene_to_trans_names = {}
    print(f"Got {len(conn_components)} connected components")
    for gene_num, component in enumerate(conn_components, 1):
        gene_id = f"reg_{gene_num}"  # give a gene some name
        # get transcripts and their ranges -> saved as component nodes names
        transcripts = set(component.nodes())
        gene_to_trans_names[gene_id] = transcripts
    return gene_to_trans_names


def save_inferred(path, gene_to_trans):
    """Save isoforms inferred from data.

    Use the same format: two tab-separated columns.
    First column: gene id, second one: isoform id.
    """
    f = open(path, "w")
    f.write("gene_id\ttrascript_id\n")
    for g, ts in gene_to_trans.items():
        for t in ts:
            f.write(f"{g}\t{t}\n")
    f.close()


def get_isoforms_data(transcript_dct, isoforms_file):
    """Add isoforms data."""
    # isoforms -> gene mapping: will be written as a class attribute
    if isoforms_file is None:
        print("\n# Isoforms file was not provided: inferring from data")
        # check networkx version
        version = nx.__version__.split(".")
        fst = int(version[0])
        snd = int(version[1]) if len(version) > 1 else 9
        num_ = float(f"{fst}.{snd}")
        if num_ < 2.4:  # old networkx version
            # notify the user about this
            msg = (
                f"Warning! You use networkx v{version}\nSplitting components with "
                f"nx.connected_component_subgraphs(), which is deprecated.\n"
                f"Please upgrade networkx to suppress this warning\n"
            )
            sys.stderr.write(msg)
        # no isoforms file provided: need to infer this
        exons_list, exon_id_to_transcript = get_exons_data(transcript_dct)
        # two exons can intersect iff: they are on the same chrom and dir
        chr_dir_to_exons = split_exons_in_chr_dir(exons_list)
        # the main part -> get exon intersections
        # create a graph of transcripts with intersected exons
        conn_graph = intersect_exons(chr_dir_to_exons, exon_id_to_transcript)
        # split graph into connected components
        # if two transcripts are in the same component -> they belong to the same gene
        components = get_graph_components(conn_graph)
        gene_to_trans = parse_components(components)
        # save_inferred(save__to, gene_to_trans) if save__to else None
    else:
        print("\n# Reading provided isoforms file")
        # just read the isoforms file
        gene_to_trans = defaultdict(
            list
        )  # one-2-many connection: one gene can have >= 1 isoform
        f = open(isoforms_file, "r")
        l_num = 0  # just in case: to point if there is an error
        skipped, added = 0, 0
        for l_num, line in enumerate(f, 1):
            line_data = line.rstrip().split("\t")
            if len(line_data) != ISOFORMS_FILE_COLS:
                err_msg = (
                    f"Isoforms file {isoforms_file} line num {l_num} corrupted: "
                    f"Expected {ISOFORMS_FILE_COLS} lines, got {len(line_data)}\n"
                )
                raise ValueError(err_msg)
            gene_id = line_data[0]
            transcript_id = line_data[1]
            if transcript_id not in transcript_dct:
                # this transcript didn't appear in our input file, no need to save it
                skipped += 1
                continue
            added += 1
            gene_to_trans[gene_id].append(transcript_id)
        f.close()
        print(f"Got {l_num} lines of isoforms data")
        print(f"Added {added} isoforms, skipped: {skipped}")

    # write isoforms data to trans objects
    for gene_id, transcript_ids in gene_to_trans.items():
        for transcript_id in transcript_ids:
            transcript_dct[transcript_id].corresponding_gene = gene_id
    # check whether there are genes without isoforms
    trans_not_assigned = []
    if isoforms_file:
        # if isoforms file is provided: there is a chance that there is some isoform
        # that does not have a corresponding gene
        for t_obj in transcript_dct.values():
            if t_obj.corresponding_gene is None:
                print_stderr(f"Warning! No gene assigned to {t_obj.id}")
                trans_not_assigned.append(t_obj.id)
    # if isoforms inferred from data: all coding isoforms must have a gene associated
    # non-coding are skipped automatically: we use CDS to infer genes
    print("#\n")
    return gene_to_trans, trans_not_assigned  # get_isoforms_data


def invert_complement(seq):
    """Make inverted-complement sequence."""
    reverse = seq[::-1]  # put N for unknown characters:
    reverse_complement = "".join(
        [complement.get(c) if complement.get(c) else "N" for c in reverse]
    )
    return reverse_complement


def split_trans_into_chroms(transcript_dct):
    """Create a dict chromosome: transcripts."""
    ret = defaultdict(list)
    for t_id, t_obj in transcript_dct.items():
        ret[t_obj.chrom].append(t_id)
    return ret


def add_seq_data(transcript_dct, two_bit, v=False):
    """Get data about reading frame."""
    if two_bit is None:
        # in this case we cannot extract any sequences
        print("2bit file was not provided")
        print_stderr("\n!!!Warning!!!\nTwo bit file was not provided")
        print_stderr("The filter quality might be affected")
        return
    genome = TwoBitFile(two_bit)
    print_stderr("Running sequence analysis...")
    # task_size = len(transcript_dct.keys())
    num = 1
    # split into chroms: read two bit faster
    chrom_to_trans_ids = split_trans_into_chroms(transcript_dct)
    chrom_num = len(chrom_to_trans_ids.keys())
    for num, (chrom, trans_ids) in enumerate(chrom_to_trans_ids.items(), 1):
        # go chrom-by-chrom
        chrom_seq = genome[chrom]
        t_num = len(trans_ids)
        print_stderr(
            f"Processing chrom: {chrom} ({t_num} transcripts) {num}/{chrom_num}"
        ) if v else None
        for t_id in trans_ids:
            t_obj = transcript_dct[t_id]
            if t_obj.is_coding is not True:
                # if there is no coding sequence -> we cannot extract it
                continue
            cds_part = chrom_seq[t_obj.trans_start : t_obj.trans_end].upper()
            t_obj.seq_analysis(cds_part)
        # for t_id in trans_ids:
    # mail loop end
    return  # add_seq_data


def mark_longest_isoforms(gene_to_iforms, transcipt_dct, allow_tga=False):
    """For each gene, get the longest intact isoform.

    If isoform is broken but is longer than the longest intact one,
    we mark is_longest_intact_cds_is_gene as False.
    """
    print("Getting the longest isoforms...")
    for _, isoforms in gene_to_iforms.items():
        # _ -> not really necessary here gene identifier
        # all [isoforms] belong to the same gene
        # sort isoforms by length
        intact_isoforms = [x for x in isoforms if transcipt_dct[x].is_intact()]
        if len(intact_isoforms) > 0:  # choose between intact isoforms
            use__list_of_ifoms = intact_isoforms
        else:  # all isoforms are not intact, use what we have
            use__list_of_ifoms = isoforms
        iform_lengts = sorted(
            [(x, transcipt_dct[x].tot_cds_len) for x in use__list_of_ifoms],
            key=lambda x: x[1],
            reverse=True,
        )
        # reversed order, 0 contains the maximal length
        max_len = iform_lengts[0][1]
        # if there are > 1 isoform with the max length then we have multiple longest isoforms
        # mark them as such
        for t_, l_ in iform_lengts:
            if l_ == max_len:
                transcipt_dct[t_].is_longest_intact_cds_in_gene = True
            else:
                transcipt_dct[t_].is_longest_intact_cds_in_gene = False
    return  # mark_longest_isoforms


def add_principal_isoforms_data(transcript_dct, proncipal_isoforms_file):
    """Add principal isoforms data."""
    if proncipal_isoforms_file is None:
        # nothing provided
        print("No principal isoforms provided, using the longest CDS as principal")
        return
    num_of_p_isoforms = 0
    f = open(proncipal_isoforms_file, "r")
    for line in f:
        item = line.rstrip()
        if not item in transcript_dct:
            # skip this item: not in the transcript dict
            continue
        transcript_dct[item].is_principal = True
    for v in transcript_dct.values():
        if v.is_principal is None:
            v.is_principal = False
    print(f"Added {num_of_p_isoforms} principal isoforms")
    f.close()


def mark_long_enough_isoforms(transcript_dct, gene_to_isoforms, min_len_param):
    """Mark long enough isoforms.

    Long enough (by default): an isoform that is not shorter
    than 80% of the shortest principal isoform, or, if no
    principal isoforms provided, 80% of the longest (CDS)
    isoform.
    """
    _min_perc = min_len_param * 0.01
    for _, transcirpts in gene_to_isoforms.items():
        trans_objects = [transcript_dct[t_] for t_ in transcirpts]
        principal_lens = [
            t_.tot_cds_len for t_ in trans_objects if t_.is_principal is True
        ]
        min_principal_len = 0 if len(principal_lens) == 0 else min(principal_lens)
        # will crash if this gene doesn't have the longest isoform, but probably this is not possible
        biggest_len = [
            t_.tot_cds_len
            for t_ in trans_objects
            if t_.is_longest_intact_cds_in_gene is True
        ][0]
        border = min_principal_len if min_principal_len > 0 else biggest_len
        threshold = _min_perc * border
        for t_ in transcirpts:
            if transcript_dct[t_].tot_cds_len < threshold:
                transcript_dct[t_].is_long_enough = False
            else:
                transcript_dct[t_].is_long_enough = True


def mark_non_unique_isoforms(transcript_dct, gene_to_isoforms):
    """In each gene mark non-unique CDS transcripts."""
    for _, transcripts in gene_to_isoforms.items():
        # go gene-by-gene
        trans_objects = [transcript_dct[t_] for t_ in transcripts]
        trans_cds_hashes = [t.get_cds_hash() for t in trans_objects]
        if len(set(trans_cds_hashes)) == len(trans_objects):
            # all CDS of transcripts of the gene are unique
            continue
        # there are non-unique hashes
        uniq_cds_to_ids = defaultdict(list)
        for h, t in zip(trans_cds_hashes, trans_objects):
            uniq_cds_to_ids[h].append(t.id)
        to_resolve = [k for k, v in uniq_cds_to_ids.items() if len(v) > 1]
        for h in to_resolve:
            # non-unique CDS - by - non-unique CDS
            # pick one random
            t_ids = uniq_cds_to_ids[h]
            loci_sizes = [transcript_dct[t].gene_locus_size for t in t_ids]
            # there is a chance that longest UTR isoform is also NMD
            nmd_feat = [transcript_dct[t].is_nmd__no_seq for t in t_ids]
            # so we select longest-UTR non-NMD isoform
            # but if ALL are NMD: then just select longest UTR
            size_sorted = sorted(
                zip(t_ids, loci_sizes, nmd_feat), key=lambda x: x[1], reverse=True
            )
            if all(nmd_feat):
                # ALL are NMD: ok, select the longest and that's it
                picked = size_sorted[0][0]  # [0]: elem with longest locus
            else:
                # list is loci-length sorted: no need to keep transcripts
                # select non-NMD
                # [0] -> transcript ID
                # [2] -> NMD or not
                non_nmd_sorted = [x[0] for x in size_sorted if x[2] is False]
                picked = non_nmd_sorted[0]
            # 0 -> from pair (transcript_id, locus lenght)
            for t_id in t_ids:
                if t_id == picked:
                    continue
                transcript_dct[t_id].non_uniq_replaced_with = picked
            # for h in to_resolve:
    return


def write_seq_stats(transcript_dct):
    """Get some transcipt statistics."""
    not_start = len(
        [t_ for t_ in transcript_dct.values() if t_.starts_with_atg is False]
    )
    not_end = len([t_ for t_ in transcript_dct.values() if t_.ends_with_stop is False])
    has_inframe_stops = len(
        [t_ for t_ in transcript_dct.values() if len(t_.inframe_stop_codons) > 0]
    )
    not_div_3 = len(
        [t_ for t_ in transcript_dct.values() if t_.frame_is_correct is False]
    )
    # non_coding = len([t_ for t_ in transcript_dct.values() if t_.is_coding is False])
    has_broken_introns = len(
        [t_ for t_ in transcript_dct.values() if len(t_.broken_introns) > 0]
    )

    print(
        "### Some transcript statistics: (a single transcipt may fall into several categories)"
    )
    print(f"Not start with ATG: {not_start}")
    print(f"Not end with a stop codon: {not_end}")
    print(f"Have premature stop codons: {has_inframe_stops}")
    print(f"Have non-canonical splice sites: {has_broken_introns}")
    print(f"CDS length modulo 3 != 0: {not_div_3}")
    # print(f"Are non-coding at all: {non_coding}")
    print()  # write_seq_stats


def read_list_of_transcripts(in_file):
    """List of transcripts that must be (/not) included."""
    if in_file is None:
        return set()
    with open(in_file, "r") as f:
        trans = set(l.rstrip() for l in f)
    return trans


def filter_step_1(
    transcipt_dct,
    min_intron_size,
    must_incl,
    must_del,
    soft_del,
    keep_nmd,
    keep_sec,
    keep_non_div_3,
    keep_incomplete,
    isof_len,
    uniq,
):
    """Make decisions about some of the transcripts."""
    print("\n# Main filter in progress...")
    filter_stats = defaultdict(list)
    genes_represented = []  # to keep genes that have selected isoforms
    for transcript_obj in transcipt_dct.values():
        repr_gene = transcript_obj.corresponding_gene
        multi_exon = transcript_obj.min_intron_size > 0
        micro_introns_found = transcript_obj.min_intron_size < min_intron_size

        if transcript_obj.id in must_incl:
            # user explicitly asked to add this transcript without any questions
            transcript_obj.selected = True
            genes_represented.append(repr_gene)
            continue
        # standard filtering procedure
        if transcript_obj.id in must_del:
            # user asked to delete this gene
            transcript_obj.selected = False
            filter_stats["in_must_del_list"].append(transcript_obj.id)
            continue
        # also exclude soft delete transcripts
        if transcript_obj.id in soft_del:
            transcript_obj.selected = False
            filter_stats["in_soft_del_list"].append(transcript_obj.id)
        elif transcript_obj.is_coding is False:
            # if not-coding: no way
            transcript_obj.selected = False
            filter_stats["non_coding_transcript"].append(transcript_obj.id)
            continue
        elif transcript_obj.frame_is_correct is False and keep_non_div_3 is False:
            # no frame - no way
            transcript_obj.selected = False
            filter_stats["frame len mod 3 != 0"].append(transcript_obj.id)
            continue
        elif micro_introns_found and multi_exon:
            # there is a very short intron: remore this
            transcript_obj.selected = False
            filter_stats["has micro introns"].append(transcript_obj.id)
            continue
        elif transcript_obj.is_long_enough is False:
            # is too short
            transcript_obj.selected = False
            filter_stats[f"shorter than {isof_len}% CDS of principal isoform"].append(
                transcript_obj.id
            )
            continue
        elif transcript_obj.is_nmd__no_seq is True and keep_nmd is False:
            transcript_obj.selected = False
            filter_stats["potential NMD target"].append(transcript_obj.id)
            continue
        elif transcript_obj.non_uniq_replaced_with is not None:
            transcript_obj.selected = False
            filter_stats["non-unique CDS"].append(transcript_obj.id)
            continue
        # ok, let's ignore broken introns for now
        # elif len(transcript_obj.broken_introns) > 0:
        #     transcript_obj.selected = False
        #     continue

        # check that ATG, stops etc params are False, not None
        # if None -> it was not checked
        atg_issue = transcript_obj.starts_with_atg is False
        stop_issue = transcript_obj.ends_with_stop is False
        all_stops_tga = transcript_obj.all_stops_are_tga
        inframe_stop = len(transcript_obj.inframe_stop_codons) > 0
        if keep_sec is True and all_stops_tga is True:
            # if we asked to ignore TGA inframe codons and at the same time
            # all inframe stops are TGA then correct inframe_stop variable,
            # set it to False
            inframe_stop = False
        if atg_issue and keep_incomplete is False:
            transcript_obj.selected = False
            filter_stats["no ATG at CDS start"].append(transcript_obj.id)
            continue
        elif stop_issue and keep_incomplete is False:
            transcript_obj.selected = False
            filter_stats["no stop codon at CDS end"].append(transcript_obj.id)
            continue
        elif inframe_stop:
            transcript_obj.selected = False
            filter_stats["has premature stop codons"].append(transcript_obj.id)
            continue
        # ok, end of all filters
        transcript_obj.selected = True
        genes_represented.append(repr_gene)

    # give some stats
    added, skipped = 0, 0
    for t_ in transcipt_dct.values():
        if t_.selected is True:
            added += 1
        else:
            skipped += 1
    print(f"Passed the filter: {added} transcripts")
    print(f"Skipped: {skipped} transcripts")

    return set(genes_represented), filter_stats


def write_filter_stats(stats):
    """Just write filter stats."""
    print("\n### Filter statistics")
    print("Reason:	number of transcripts omitted because of: (first exclusion reason)")
    for reason, transcripts in stats.items():
        print(f"{reason}:\t{len(transcripts)}")
    print("####\n")  # \n
    return  # write_filter_stats


def filter_step_2(transcript_dct, skipped_genes, gene_to_isoforms, must_del):
    """Trying to avoid missing genes.

    These genes have only bad isoforms, need to find at least not terrible.
    """
    genes_passed = []
    print(f"Trying to save {len(skipped_genes)} genes")
    transcripts_passed = 0
    for gene in skipped_genes:
        transcripts = gene_to_isoforms[gene]
        # there is no sense in adding transcripts that are non-coding at all
        longest_one = [
            t_ for t_ in transcripts if transcript_dct[t_].is_longest_intact_cds_in_gene
        ]
        longest_filt = [t_ for t_ in longest_one if t_ not in must_del]
        if len(longest_filt) == 0:
            continue
        for x in longest_filt:
            transcript_dct[x].selected_to_avoid_skipping_gene = True
            transcripts_passed += 1
        genes_passed.append(gene)
    print(f"Added extra {transcripts_passed} transcripts")
    return set(genes_passed)


def save_log(log_path, transcript_dict):
    """Save dataset of transcripts."""
    if log_path is None:
        return
    f = open(log_path, "w")
    header = (
        "transcript gene tot_cds_len is_coding in_frame atg_start stop_end inframe_stops "
        "u12_introns na_introns is_longest is_principal long_enough min_intron_len in_must_incl_list "
        "all_stops_tga is_nmd non_uniq selected extra_selection".split()
    )
    f.write("\t".join(header))
    f.write("\n")
    for trans in transcript_dict.values():
        data = (
            trans.id,
            trans.corresponding_gene,
            trans.tot_cds_len,
            trans.is_coding,
            trans.frame_is_correct,
            trans.starts_with_atg,
            trans.ends_with_stop,
            trans.inframe_stop_codons,
            trans.u12_introns,
            trans.broken_introns,
            trans.is_longest_intact_cds_in_gene,
            trans.is_principal,
            trans.is_long_enough,
            trans.min_intron_size,
            trans.must_include,
            trans.all_stops_are_tga,
            trans.is_nmd__no_seq,
            trans.non_uniq_replaced_with,
            trans.selected,
            trans.selected_to_avoid_skipping_gene,
        )
        data_str = "\t".join([str(x) for x in data])
        f.write(data_str)
        f.write("\n")
    f.close()


def save_out(transcript_dct, output_file, out_fmt):
    """Write transcripts to output file."""
    saved_counter = 0
    f = open(output_file, "w")
    for transcript_obj in transcript_dct.values():
        to_be_added = (
            transcript_obj.selected or transcript_obj.selected_to_avoid_skipping_gene
        )
        if to_be_added is False:
            continue
        if out_fmt == ORIGINAL:
            # save in the same format as input
            line = transcript_obj.raw_line.rstrip()
        elif out_fmt == transcript_obj.in_fmt:
            # out_fmt matches the in_fmt
            line = transcript_obj.raw_line.rstrip()
        elif out_fmt == BED:
            line = transcript_obj.output_bed()
        elif out_fmt == GENEPRED:
            line = transcript_obj.output_genepred()
        else:
            line = ""
            raise NotImplementedError(f"{out_fmt} is not yet supported for ouput")
        f.write(f"{line}\n")
        saved_counter += 1
    f.close()
    print(f"\n#Done, saved {saved_counter} transcripts")


def write_warnings(t_without_gene, skipped_genes, added_skipped_genes, wd):
    """Write useful warning message at the end."""
    if len(t_without_gene) > 0:
        # there are some transcripts without associated gene
        # this is because of lack of data in the isoforms file provided by user
        print_stderr("\n!!!WARNING!!!")
        print_stderr("The provided isoforms file is incomplete")
        print_stderr("The following transcripts don't have an associated gene:")
        for t in t_without_gene:
            print_stderr(t)
        print_stderr(
            "Note: using the results as input for TOGA will result in error exit"
        )
        print_stderr("because TOGA requires each transcript to have an associated gene")
        print_stderr("######")
    skipped_genes_num = len(skipped_genes)
    if skipped_genes_num > 0 and added_skipped_genes is False:
        # we skipped some genes completely: need to notify user
        print_stderr("\n!!!WARNING!!!")
        print_stderr(f"{skipped_genes_num} genes were skipped")
        write_to_console = skipped_genes_num < MAX_NUM_SKIPPED_GENES_TO_CONSOLE
        # --> If >50, pls write them to a 'excluded.genes.list.txt' file instead of printing them all out. ]
        if write_to_console:
            for gene in skipped_genes:
                # Well, there is an issue: transcripts can be skipped due to various reasons
                # there is no "a single reason" why some gene is omitted
                # instead, it's like transcript 1 skipped because reason a b and c
                # transcript 2 because of reasons b and d, et cetera
                print_stderr(gene)
        else:
            # write to a file
            time_ = dt.now().strftime("%Y_%m_%d_at_%H_%M_run")  # to avoid collisions
            filename = f"{time_}_excluded_genes.txt"
            path = os.path.join(wd, filename)
            print_stderr(f"List of the excluded genes saved to {path}")
            f = open(path, "w")
            for gene in skipped_genes:
                f.write(f"{gene}\n")
            f.close()

        print_stderr(
            "All isoforms of these genes had issues critical enough to filter them out"
        )
        print_stderr("If you need to include these genes, call the script with -k flag")
        print_stderr(
            "Then the script will include the longest CDS isoform of these genes"
        )
        print_stderr("Pls also consider adding transcripts with premature TAG codons")
        print_stderr("These transcripts may encode selenocysteine, to do so add")
        print_stderr("--sec flag to the command call")
        print_stderr("####")
    elif len(skipped_genes) > 0 and added_skipped_genes is True:
        # we would've skipped some genes but we forced the scipt to add them
        pass
    # TODO: check whether anything else is needed


def parse_args():
    """Read args, check."""
    app = argparse.ArgumentParser(description="Filter transcripts for TOGA input.")
    app.add_argument(
        "input", help="Annotation file to be filtered (bed or genePredEx)."
    )
    app.add_argument("output", help="Path to the filtered annotation file.")
    app.add_argument(
        "--in_fmt",
        default="bed",
        help="Format of input file: IN_FMT bed -> bed-12 format "
        "IN_FMT genePred -> genepred format (15 columns), "
        "default: bed",
    )
    app.add_argument("--two_bit", default=None, help="Path to the genome 2bit file")
    app.add_argument(
        "--out_fmt",
        default=ORIGINAL,
        help="Output format, by default save in the same " "format as input",
    )
    app.add_argument(
        "--isoforms",
        "-i",
        default=None,
        help="Annotation isoforms data. tab-delim file containing 2 columns, where "
        "1st column: gene_id, 2nd column: transcript_id. If not provided, "
        "the script will infer isoforms from data",
    )
    app.add_argument(
        "--min_isoform_length",
        type=int,
        default=80,
        help="Isoforms that represent at least <minIsoformLength> percentge of "
        "the longest transcript will be kept. Values from 1 to 100 (default) "
        "--minIsoformLength 80",
    )
    app.add_argument(
        "--principal_isoforms",
        "--pi",
        default=None,
        help="If provided, use as must include isoforms.",
    )
    app.add_argument(
        "--include_sec_transcripts",
        "--sec",
        dest="include_sec_transcripts",
        action="store_true",
        help="Do not exclude transcripts with premature TGA codons, "
        "which can be selenocysteine codons. Transcripts with other "
        " stop codons (TAA, TAG) will still be filtered out. ",
    )
    app.add_argument(
        "--keep_nmd",
        default=False,
        action="store_true",
        dest="keep_nmd",
        help="keep targets of Nonsense-Mediated mRNA Decay (not recommended).",
    )
    app.add_argument(
        "--min_intron_size",
        "--mi",
        type=int,
        default=20,
        help="Filter out transcripts with introns shorter than <minIntronSize>. "
        "Default value is 20",
    )
    app.add_argument(
        "--uniq",
        "-u",
        action="store_true",
        dest="uniq",
        help="Keep transcripts with unique CDS only",
    )
    app.add_argument(
        "--keep_broken_transcripts",
        "-k",
        default=False,
        action="store_true",
        dest="keep_broken_transcripts",
        help="If set, keep broken transcripts (those which would be filtered out "
        "by NMD, intron size, incorrect reading frame, etc.) if all transcripts "
        "of a gene would be filtered out. This avoids skipping a gene entirely.",
    )
    app.add_argument(
        "--keep_incomplete_transcripts",
        "--kit",
        default=False,
        action="store_true",
        dest="keep_incomplete_transcripts",
        help="If set, keep transcripts with CDS length modulo 3 != 0"
    )
    app.add_argument(
        "--filter_dataframe",
        "--fdf",
        default=None,
        help="Save features of each transcript used for classification  (strictly recommended)",
    )
    app.add_argument(
        "--deleted_genes_list",
        "--dgf",
        default=None,
        help="Output entirely deleted genes",
    )
    app.add_argument(
        "--must_include_list",
        default=None,
        help="A text file containing a list of transcripts " "that must be kept",
    )
    app.add_argument(
        "--allow_incomplete_ends",
        action="store_true",
        dest="allow_incomplete_ends",
        help="Keep transcripts lacking complete 5' or 3' ends: not starting with ATG and/or "
        "not ending with a stop codon",
    )
    app.add_argument(
        "--must_delete_list",
        default=None,
        help="List of transcripts that must be deleted",
    )
    app.add_argument(
        "--sort_del_list",
        "--sd",
        default=None,
        help=(
            "List of transcripts to be deleted, the difference with must_delete list: in case "
            " --keep_broken_transcripts flag is on, these transcrips will be recovered"
        ),
    )
    app.add_argument(
        "--save_isoforms",
        default=None,
        help="If isoforms were inferred from data by "
        "same-strand CDS overlap save them into "
        "the given output file",
    )
    app.add_argument(
        "--verbose", "-v", action="store_true", dest="verbose", help="Enable verbosity"
    )

    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    # args sanity checsks
    if (
        args.in_fmt.lower() not in SUPPORTED_IN_FMT
    ):  # check that input format is supported
        raise ValueError(
            f"Error! --in_fmt expected to be in {SUPPORTED_IN_FMT}, got {args.in_fmt}"
        )
    elif (
        args.out_fmt.lower() not in SUPPORTED_OUT_FMT
    ):  # also check that output format is supported
        raise ValueError(
            f"Error! --out_fmt expected to be in {SUPPORTED_OUT_FMT}, got {args.out_fmt}"
        )
    elif (
        args.min_isoform_length < 0 or args.min_isoform_length > 100
    ):  # percents, in 0..100
        raise ValueError(
            f"min_isoform length arument must be in [0..100], got {args.min_isoform_length}"
        )
    elif not os.path.isfile(args.input):  # check that input file exists
        raise FileNotFoundError(
            f"Is not a file: {args.input} (plt check input argument)"
        )
    # if any of these arguments set: also check that files exist
    elif args.isoforms and not os.path.isfile(args.isoforms):
        raise FileNotFoundError(
            f"Is not a file: {args.isoforms} (pls check --isoforms argument)"
        )
    elif args.must_include_list and not os.path.isfile(args.must_include_list):
        raise FileNotFoundError(
            f"Is not a file: {args.must_include_list} (pls check --must_include_list argument)"
        )
    elif args.must_delete_list and not os.path.isfile(args.must_delete_list):
        raise FileNotFoundError(
            f"Is not a file: {args.must_delete_list} (pls check --must_delete_list argument)"
        )
    elif args.principal_isoforms and not os.path.isfile(args.principal_isoforms):
        raise FileNotFoundError(
            f"Is not a file: {args.principal_isoforms} (pls check --principal_isoforms argument)"
        )

    return args


def save_deleted_genes_data(deleted_genes_list, skipped_genes, gene_to_isoforms):
    """In case some genes are entirely deleted list them (is user asks for)."""
    if deleted_genes_list is None:
        return
    f = open(deleted_genes_list, "w")
    f.write("GENE\tTRANSCRIPT\n")
    for gene in skipped_genes:
        isoforms = gene_to_isoforms.get(gene)
        for isoform in isoforms:
            f.write(f"{gene}\t{isoform}\n")
    f.close()


def save_isoforms(save_to, gene_to_isoforms):
    """Save isoforms to a given file.
    
    Save only those transcripts that passed the filter."""
    f = open(save_to, "w")
    f.write("GENE\tTRANSCRIPT\n")
    for gene, isoforms in gene_to_isoforms.items():
        for isoform in isoforms:
            f.write(f"{gene}\t{isoform}\n")
    f.close()


def select_genes_and_transcripts_to_save(transcript_dct):
    gene_to_isoforms_final = defaultdict(list)
    for t_id, transcript_obj in transcript_dct.items():
        to_be_added = (
        transcript_obj.selected or transcript_obj.selected_to_avoid_skipping_gene
        )
        if to_be_added is False:
            continue
        gene_id = transcript_obj.corresponding_gene
        if gene_id is None:
            gene_id = t_id
        gene_to_isoforms_final[gene_id].append(t_id)
    return gene_to_isoforms_final


def main():
    """Entry point."""
    args = parse_args()
    # read lists of must_include and must_delete genes (if provided)
    must_incl = read_list_of_transcripts(args.must_include_list)
    must_del = read_list_of_transcripts(args.must_delete_list)
    soft_del = read_list_of_transcripts(args.sort_del_list)

    i_d_intersect = must_incl.intersection(must_del)
    # also check that these lists do not intersect
    if len(i_d_intersect) > 0:
        err_msg = (
            "Please check your must delete and must include transcript sets "
            "because they do intersect; transcripts that appear in both lists:"
        )
        print_stderr(err_msg)
        for elem in i_d_intersect:
            print_stderr(elem)
        raise ValueError("Same gene in must include and must delete list")
    
    must_and_soft_del = must_del.intersection(soft_del)
    if len(must_and_soft_del) > 0:
        err_msg = (
            f"Must_ and soft_ delete sets have {len(must_and_soft_del)} transcripts "
            f"in common! Please exclude them to avoid ambiguity"
        )
    # read annotation file and create transcripts dictionary
    transcript_dct = read_transcripts(args.input, args.in_fmt)
    # read isoforms data if provided, otherwise infer isoforms from data
    gene_to_isoforms, t_without_gene = get_isoforms_data(
        transcript_dct, args.isoforms
    )
    # add sequences data (start/stop codons, etc) if 2bit is provided
    add_seq_data(transcript_dct, args.two_bit, v=args.verbose)
    if args.two_bit:
        write_seq_stats(transcript_dct)
    # get longest isoforms per gene + principal (if provided)
    mark_longest_isoforms(
        gene_to_isoforms, transcript_dct, allow_tga=args.include_sec_transcripts
    )
    add_principal_isoforms_data(transcript_dct, args.principal_isoforms)
    mark_long_enough_isoforms(transcript_dct, gene_to_isoforms, args.min_isoform_length)
    mark_non_unique_isoforms(transcript_dct, gene_to_isoforms)

    genes_represented, fstat = filter_step_1(
        transcript_dct,
        args.min_intron_size,
        must_incl,
        must_del,
        soft_del,
        args.keep_nmd,
        args.include_sec_transcripts,
        args.keep_incomplete_transcripts,
        args.allow_incomplete_ends,
        args.min_isoform_length,
        args.uniq,
    )
    write_filter_stats(fstat)
    all_genes = set(gene_to_isoforms.keys())
    print(f"Total number of genes: {len(all_genes)}")
    print(f"Added isoforms of {len(genes_represented)} genes")
    skipped_genes = all_genes.difference(genes_represented)
    print(f"Genes not represented: {len(skipped_genes)}")
    if len(skipped_genes) > 0 and args.keep_broken_transcripts:
        # we like to avoid skipping genes
        # in this case, we try to include as good transcripts as possible
        genes_represented_2 = filter_step_2(
            transcript_dct, skipped_genes, gene_to_isoforms, must_del
        )
        print(f"Included {len(genes_represented_2)} genes in addition")
        skipped_genes = skipped_genes.difference(genes_represented_2)
    # save the filtered items
    save_log(args.filter_dataframe, transcript_dct)
    save_deleted_genes_data(args.deleted_genes_list, skipped_genes, gene_to_isoforms)
    save_out(transcript_dct, args.output, args.out_fmt)
    wd = os.path.dirname(
        args.output
    )  # to save all necessary warnings in the same dir with input
    if args.save_isoforms:
        kept_gene_to_isoforms = select_genes_and_transcripts_to_save(transcript_dct)
        save_isoforms(args.save_isoforms, kept_gene_to_isoforms)
    write_warnings(t_without_gene, skipped_genes, args.keep_broken_transcripts, wd)


if __name__ == "__main__":
    main()

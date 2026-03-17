from datetime import datetime
from lib.deletion import DeletionSet, read_all_line_subs
import subprocess
from subprocess import CalledProcessError
import logging

HEADERS = f"""##fileformat=VCFv4.2
##ALT=<ID=DEL,Description="Deletion">
##fileDate={datetime.now().strftime("%Y%m%d")}
##source=script generation
##reference=Pisum_sativum-JI2822-JIC_v1.3.fasta
##contig=<ID=chr1,length=480540775>
##contig=<ID=chr2,length=507322098>
##contig=<ID=chr3,length=538651075>
##contig=<ID=chr4,length=511007960>
##contig=<ID=chr5,length=647395306>
##contig=<ID=chr6,length=535013209>
##contig=<ID=chr7,length=563422397>
##commandline="uv run main.py"
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of deletion.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of deletion.">
##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""

logger = logging.getLogger(__name__)

CHROMOSOME_NAME_MAPPING = {f"chr{i + 1}": f"OZ0754{28 + i}.1" for i in range(7)}


def write_vcf_file(merged_deletion_sets: list[DeletionSet]) -> None:
    logger.info("Creating VCF file...")
    with open("output.vcf", "w") as f:
        f.write(HEADERS + "\n")
        all_line_subs = read_all_line_subs()
        f.write(
            "#"
            + "\t".join(
                ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                + list(all_line_subs)
            )
            + "\n"
        )

        # mapping of deletion chromosome, start and end point a mapping of lines it is present in to whetehr it is homoyzgous
        lines: dict[tuple[str, int, int], dict[str, bool]] = {}
        for merged_deletion_set in merged_deletion_sets:
            line_sub = merged_deletion_set.line_sub
            chromosome = merged_deletion_set.chromosome
            for deletion in merged_deletion_set.deletions:
                key = (chromosome, deletion.start, deletion.end)
                if key in lines:
                    if line_sub in lines[key]:
                        raise ValueError("How is there a duplicate")
                    lines[key][line_sub] = deletion.homozygous
                else:
                    lines[key] = {line_sub: deletion.homozygous}

        # sort by chromsome then by start point
        lines = dict(
            sorted(
                lines.items(), key=lambda x: (int(x[0][0].removeprefix("chr")), x[0][1])
            )
        )
        for (chromosome, start, end), line_sub_mapping in lines.items():
            chromosome_silly_name = CHROMOSOME_NAME_MAPPING[chromosome]
            length = end - start
            nucleotide_query = subprocess.run(
                [
                    "samtools",
                    "faidx",
                    "Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa",
                    f"{chromosome_silly_name}:{start}-{start}",
                ],
                capture_output=True,
                text=True,
            )
            try:
                nucleotide_query.check_returncode()
            except CalledProcessError:
                logger.error(nucleotide_query.stderr)
            nucleotide = nucleotide_query.stdout.splitlines()[1]
            values = [
                chromosome,
                start,
                ".",
                nucleotide,
                "<DEL>",
                ".",
                ".",
                ";".join((f"SVLEN={-length}", "SVCLAIM=D", f"END={end}")),
                "GT",
            ]
            for line_sub_iter in all_line_subs:
                if line_sub_iter in line_sub_mapping:
                    values.append("1/1" if line_sub_mapping[line_sub_iter] else "0/1")
                else:
                    values.append("0/0")
            f.write("\t".join(map(str, values)) + "\n")

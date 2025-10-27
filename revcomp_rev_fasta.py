#!/usr/bin/env python3
"""
reverse_fasta.py

对FASTA文件中的参考序列执行反向互补或仅反向操作。
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    """返回DNA序列的反向互补序列"""
    complement = str.maketrans("ATGCatgcNn", "TACGtacgNn")
    return seq.translate(complement)[::-1]

def reverse_fasta(input_fasta, output_fasta, mode="reverse_complement"):
    """对FASTA文件执行反向互补或反向操作"""
    records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq)
        if mode == "reverse_complement":
            new_seq = reverse_complement(seq)
            suffix = "_revcomp"
        elif mode == "reverse":
            new_seq = seq[::-1]
            suffix = "_rev"
        else:
            raise ValueError("mode 必须是 'reverse' 或 'reverse_complement'")

        record.id = record.id + suffix
        record.description = ""
        record.seq = Seq(new_seq)   # ✅ 关键修改
        records.append(record)

    SeqIO.write(records, output_fasta, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="对FASTA序列进行反向互补或反向")
    parser.add_argument("input_fasta", help="输入fasta文件路径")
    parser.add_argument("output_fasta", help="输出fasta文件路径")
    parser.add_argument("--mode", choices=["reverse", "reverse_complement"], default="reverse_complement",
                        help="选择操作模式：reverse 或 reverse_complement (默认: reverse_complement)")
    args = parser.parse_args()

    reverse_fasta(args.input_fasta, args.output_fasta, args.mode)
    print(f"✅ 已生成 {args.output_fasta}，模式：{args.mode}")

from Bio import SeqIO

def trim_seq(in_fname, out_fname, left_clip, right_clip, min_len):

        seqs_trimmed = []
        for seq in SeqIO.parse(in_fname, 'fasta'):
            left = left_clip
            right = len(seq) - right_clip
            if right <= left:
                continue
            elif (right - left) < min_len:
                continue
            else:
                seq = seq[left:right]
                seqs_trimmed.append(seq)
        SeqIO.write(seqs_trimmed, out_fname, 'fasta')

if __name__ == '__main__':
    in_fname = 'fastqjoin.join.fasta'
    out_fname = 'fastqjoin.join.fasta'
    left_clip = 21
    right_clip = 20
    min_len = 100
    trim_seq(in_fname, out_fname, left_clip, right_clip, min_len)



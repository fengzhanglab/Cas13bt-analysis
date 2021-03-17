import Bio.SeqIO
import os
import pickle

spacers = []
with open('spacers.csv', 'r') as f:
    for line in f:
        if line.strip() == '':
            continue
        spacers.append(line.strip())
spacers = set(spacers)

l_flank = 'TAGCTCAGTCCTAGGTATAATGCTAGC'

samples = ['Bt1_input_S7_R1_001.fastq']

for sample in samples:
    total_count = 0
    count = 0
    spacer_dict = {u : 0 for u in spacers}
    with open(sample, 'r') as f:
        # Iterate through fastq file and extract records matching the template seed
        records = Bio.SeqIO.parse(f, 'fastq')
        for record in records:
            if total_count % 1000000 == 0:
                print(count, total_count, float(count) / (total_count + 1e-5))
                print(len([1 for u,v in spacer_dict.items() if v >= count*0.000001 and v >= 20]))
                print(max(spacer_dict.values()))

            seq = str(record.seq)
            quality = record.letter_annotations["phred_quality"]
            min_quality = min(quality)
            total_count += 1
            # Throw out low quality sequences
            if min_quality < 10:
                continue
            l_idx = seq.find(l_flank)
            if l_idx == -1:
                continue
            q_sub = quality[l_idx+len(l_flank):l_idx+len(l_flank)+spacer_len]
            if len(q_sub) < spacer_len:
                continue
            q = min(q_sub)
            # Throw out low quality sequences in the PAM/PFS region
            if q < 20:
                continue
            spacer = seq[l_idx+len(l_flank):l_idx+len(l_flank)+spacer_len]
                                                                                                                                                                                                                         
            if not spacer in spacer_dict:
                continue
            else:
                spacer_dict[spacer] += 1
            count += 1

    with open(os.path.splitext(sample)[0]+'spacers.p', 'wb') as f:
        pickle.dump(spacer_dict, f)








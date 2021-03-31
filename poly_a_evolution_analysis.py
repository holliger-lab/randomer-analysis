import sys, os, re
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
import seaborn as sns
from pylab import *
from sklearn.decomposition import PCA
# from sklearn.manifold import TSNE

# from openTSNE import TSNE
from tsnecuda import TSNE
import random

plt.style.use('ggplot')
sns.set_palette("husl")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def readFasta(fastaFile):
    fh = open(fastaFile, 'r')
    for line in fh:
        header = ""
        seq = ""
        if line[0] == '>':
            header = line.rstrip()[1:]
            seq = fh.readline().rstrip()
        yield [header, seq]
    fh.close()


def calcShannonIndex(uniqueSet):
    n_unique = len(uniqueSet)
    pool_size = np.sum([c for c in uniqueSet.values()])
    avg_seq_len = np.mean([len(c) for c in uniqueSet.keys()])

    H = 0.0

    for unique in uniqueSet.items():
        pi = float(int(unique[1])/pool_size)
        H += pi*np.log2(pi)

    print("Unique sequences: %d\tTotal sequences: %d\nShannon index: %f\tMean sequence length: %f\n" % (n_unique, pool_size, -H, avg_seq_len))

    return -H

def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def plot_entropy_over_rounds(fileList, bc_1, bc_2):
    '''Plots the information entropy of multiple rounds of sequences.'''

    for file in fileList:
        unique_sequences = {}
        print(file)

        pat_start = r"(" + bc_1 + ")"
        pat_start = re.compile(pat_start)

        pat_end = r"(" + bc_2 + ")"
        pat_end = re.compile(pat_end)

        for h,seq in readFasta(file):

            search_start = pat_start.search(seq)
            search_end = pat_end.search(seq)

            if search_start and search_end:
                seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
                if len(seq_trimmed) == 39:
                    if seq_trimmed in unique_sequences:
                        unique_sequences[seq_trimmed] += 1
                    else:
                        unique_sequences[seq_trimmed] = 1

        negH = calcShannonIndex(unique_sequences)



def plot_hamming_distance_over_rounds(fileList, bc_1, bc_2, ref):

    n_col = 5
    n_row = 2
    plt.figure(figsize=(15, 6))

    for i, file in enumerate(fileList):
        unique_sequences = {}
        print(file)

        pat_start = r"(" + bc_1 + ")"
        pat_start = re.compile(pat_start)

        pat_end = r"(" + bc_2 + ")"
        pat_end = re.compile(pat_end)

        for h,seq in readFasta(file):

            search_start = pat_start.search(seq)
            search_end = pat_end.search(seq)

            if search_start and search_end:
                seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
                if len(seq_trimmed) == len(ref):
                    if seq_trimmed in unique_sequences:
                        unique_sequences[seq_trimmed] += 1
                    else:
                        unique_sequences[seq_trimmed] = 1

        hd_list = []
        hd_dist = np.zeros(len(ref))
        for uniqe_seq in unique_sequences:
            hd = hamming_distance(uniqe_seq, ref)
            hd_list.append(hd)
            hd_dist[hd] += 1

        print("Mean HD: %f" % (np.mean(hd_list)))
        print("HD distribution: ", hd_dist)

        plt.subplot(n_row, n_col, i+1)
        plt.title("A1 R%d" % (i))
        plt.xlabel("HD from 359")
        plt.ylabel("Count")
        plt.bar(np.arange(hd_dist.shape[0]), hd_dist)

    plt.tight_layout()
    plt.show(block=True)
    # plt.savefig('randomer_fig4_hd_from_A1_359.pdf', dpi=300, bbox_inches='tight')


def gen_rand_seqs(seq_len, n_seq):
    nt_arr = ['A', 'T', 'C', 'G']

    rand_seqs = []

    for i in range(n_seq):
        rand_dna = ""

        for i in range(seq_len):
            rand_dna += nt_arr[random.randint(0, 3)]

        rand_seqs.append(rand_dna)

    return rand_seqs




def plot_hamming_distance_over_rounds_stacked(fileList, bc_1, bc_2, ref, labels=None, add_rand=False, write_csv=False, output_file_name="default"):

    if write_csv:
        fh = open(output_file_name + ".csv", 'w')

        header = "Experiment"
        for i in range(len(ref)):
            header += ",%d" % (i)

        fh.write(header+'\n')

    n_col = 1
    n_row = len(fileList)
    if add_rand:
        n_row += 1

    plt.figure(figsize=(4, len(fileList)))
    c = 0
    for i, file in enumerate(fileList):
        unique_sequences = {}
        print(file)

        pat_start = r"(" + bc_1 + ")"
        pat_start = re.compile(pat_start)

        pat_end = r"(" + bc_2 + ")"
        pat_end = re.compile(pat_end)

        for h,seq in readFasta(file):

            search_start = pat_start.search(seq)
            search_end = pat_end.search(seq)

            if search_start and search_end:
                seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
                if len(seq_trimmed) == len(ref):
                    if seq_trimmed in unique_sequences:
                        unique_sequences[seq_trimmed] += 1
                    else:
                        unique_sequences[seq_trimmed] = 1

        hd_list = []
        hd_dist = np.zeros(len(ref))
        for uniqe_seq in unique_sequences:
            hd = hamming_distance(uniqe_seq, ref)
            hd_list.append(hd)
            hd_dist[hd] += 1

            # if hd <= 10:
            #     print(">%s_%d" %(file, c))
            #     print(uniqe_seq)

            c += 1

        print("Mean HD: %f" % (np.mean(hd_list)))
        print("HD distribution: ", hd_dist)

        if write_csv:
            csv_line = "%s," % (file)
            csv_line += ",".join([str(x) for x in hd_dist])
            fh.write(csv_line+"\n")


        plt.subplot(n_row, n_col, i+1)
        # plt.title("A1 R%d" % (i))
        # plt.xlabel("HD from 359")
        # plt.ylabel("Count")
        if labels != None:
            plt.bar(np.arange(hd_dist.shape[0]), hd_dist, label="%s" % (labels[i]))
        else:
            plt.bar(np.arange(hd_dist.shape[0]), hd_dist, label="R%d" % (i))

        plt.yticks([])

        if add_rand == True:
            plt.xticks([0, 10, 20, 30, 40])
            plt.tick_params(axis='x', labelsize=0, length = 3, labelbottom=False)

        elif add_rand == False:
            if i != len(fileList)-1:
                plt.xticks([0, 10, 20, 30, 40])
                plt.tick_params(axis='x', labelsize=0, length = 3, labelbottom=False)
            else:
                plt.xticks([0, 10, 20, 30, 40])
                plt.xlabel("HD from seqT1")

        plt.legend()


    if add_rand: ## Random sequences.
        print("Random")
        rand_seq = gen_rand_seqs(len(ref), 10000) # Generate 10k random seqs.
        hd_list = []
        hd_dist = np.zeros(len(ref))
        for seq in rand_seq:
            hd = hamming_distance(seq, ref)
            hd_list.append(hd)
            hd_dist[hd] += 1

        print("Mean HD: %f" % (np.mean(hd_list)))
        print("HD distribution: ", hd_dist)

        if write_csv:
            csv_line = "Random,"
            csv_line += ",".join([str(x) for x in hd_dist])
            fh.write(csv_line+"\n")

        plt.subplot(n_row, n_col, n_row)
        plt.bar(np.arange(hd_dist.shape[0]), hd_dist, label="Random")
        plt.yticks([])

        plt.legend()

        plt.xticks([0, 10, 20, 30, 40])
        plt.xlabel("HD from seqT1")

    plt.tight_layout()
    # plt.show(block=True)
    plt.savefig(output_file_name + ".pdf", dpi=300, bbox_inches='tight')




nt_arr = ['A', 'T', 'C', 'G']

def generate_random_rna(nt_len):
    rand_rna = ""

    for i in range(nt_len):
        rand_rna += nt_arr[random.randint(0, 3)]

    return rand_rna

# encode_dict = {
#     'A': 0,
#     'T': 1,
#     'C': 2,
#     'G': 3,
#     'N': 4
# }

encode_dict = {
    'A': [1,0,0],
    'T': [1,1,0],
    'C': [1,1,0],
    'G': [0,0,1],
    ' ': [0,0,0]
}

def encode(sequence):
    encoded = np.zeros((45, 3), dtype=np.int)
    for i, base in enumerate(sequence):
        encoded[i] = encode_dict[base]
    return encoded.flatten()

def figTwoManifold(fileList, bc_1, bc_2, output_file_name="default", labels=None, gen_random=True):

    X = []
    X_indices = []
    j = 0
    # labels = []

    for i, file in enumerate(fileList):
        unique_sequences = {}
        print(file)
        # labels.append(i)

        pat_start = r"(" + bc_1 + ")"
        pat_start = re.compile(pat_start)

        pat_end = r"(" + bc_2[:-3] + ")"
        pat_end = re.compile(pat_end)

        tmp_indices = [j,j]

        c = 0
        for h,seq in readFasta(file):

            if 'N' in seq:
                continue

            search_start = pat_start.search(seq)
            search_end = pat_end.search(seq)

            if search_start and search_end:
                seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
                # seq_full = bc_1 + seq_trimmed + bc_2
                seq_full = seq_trimmed
                if len(seq_full) > 45:
                    continue
                    # enc = encode(seq_trimmed)
                    # x_gp.append(enc)
                    # y_gp.append(i)

                # print(seq_trimmed)
                if seq_full in unique_sequences:
                    unique_sequences[seq_full] += 1
                else:
                    unique_sequences[seq_full] = 1
                    c+=1

            # if c > 100000:
            #     break

        for seq in unique_sequences:
            enc = encode(seq)
            X.append(enc)
            j+=1

        tmp_indices[1] = j
        X_indices.append(tmp_indices)

        print(len(unique_sequences))


    if gen_random:
        print("Random")
        tmp_indices = [j,j]

        for k in range(2000000):
            seq = generate_random_rna(random.randint(32, 42))

            enc = encode(seq)
            X.append(enc)
            j+=1

        tmp_indices[1] = j
        X_indices.append(tmp_indices)
        print(k)


    ## Add specials.
    print("Specials")
    tmp_indices = [j,j]

    wt_seq = encode("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    X.append(wt_seq)
    j+=1

    selected_seq = encode("GAAGGAAGAAAATGCAGAAAAAAAGAAAAAAATGTCTGG")
    X.append(selected_seq)
    j+=1


    for h, seq in readFasta("T5_R8_hits_strict.fasta"):
        search_start = pat_start.search(seq)
        search_end = pat_end.search(seq)

        if search_start and search_end:
            seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
            X.append(encode(seq_trimmed))
            j+=1
            print(seq_trimmed)

    tmp_indices[1] = j
    X_indices.append(tmp_indices)



    # print(X_indices[-1][0]+1)

    # exit()




    X = np.array(X)
    print(X.shape)
    print(X_indices)

    # exit()

    plt.figure(figsize=(12, 6))
    ## Look at manifolds.

    # pca = PCA(n_components=2).fit(X)
    # pca_all = TSNE(n_components=2).fit_transform(X)
    # pca = TruncatedSVD(n_components=2, n_iter=7, random_state=42).fit(X)
    # pca_all = pca.transform(X)


    # tsne = TSNE(
    # perplexity=50,
    # metric="cosine",
    # initialization="pca",
    # n_jobs=16,
    # random_state=42,
    # verbose=True
    # )
    # pca_fit = tsne.fit(X)
    # pca = pca_fit.transform(X)

    tsne = TSNE(
    perplexity=50,
    metric="euclidean",
    random_seed=42,
    verbose=True
    )
    pca = tsne.fit_transform(X)

    # pca = umap.UMAP(n_neighbors=30,
    # n_components=2,
    # random_state=42)
    # pca.fit(X)

    n_col = 5
    n_row = 2

    cmap = sns.color_palette("tab10", 8)

    for i in range(len(X_indices)-1):
        idx_grp = X_indices[i]
        plt.subplot(n_row, n_col, i+1)
        # pca_all = tsne_fit.transform(x_gp)
        # pca_all = tsne.fit_transform(x_gp, pca)

        print("Plotting indices: %d:%d" %(idx_grp[0], idx_grp[1]))

        if labels == None:
            plt.title("R%d" % (i))
        else:
            plt.title(labels[i])

        # plt.scatter(pca_all[:, 0], pca_all[:, 1], s=3, c="slategrey")
        sns.kdeplot(x=pca[idx_grp[0]:idx_grp[1], 0], y=pca[idx_grp[0]:idx_grp[1], 1], cmap="viridis", fill=True)

        # print(encode("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"))
        # wt = pca.transform(encode(bc_1 + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + bc_2).reshape(1,-1))
        # selected = pca.transform(encode(bc_1 + "GAAGGAAGAAAATGCAGAAAAAAAGAAAAAAATGTCTGG" + bc_2).reshape(1,-1))
        # wt = tsne.fit_transform(encode("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").reshape(1,-1))
        # selected = tsne.fit_transform(encode("GAAGGAAGAAAATGCAGAAAAAAAGAAAAAAATGTCTGG").reshape(1,-1))


        wt_index = X_indices[-1][0]
        selected_index = X_indices[-1][0]+1
        plt.scatter(pca[wt_index, 0], pca[wt_index, 1], s=10, c='red')
        # plt.scatter(pca[selected_index, 0], pca[selected_index, 1], s=10, c='red', marker='D')

        # motif_matches_indices = [X_indices[-1][0]+2, X_indices[-1][1]-1]
        # print(motif_matches_indices)
        # plt.scatter(pca[motif_matches_indices[0]:motif_matches_indices[1], 0], pca[motif_matches_indices[0]:motif_matches_indices[1], 1], s=10, c='pink', marker='P')


        plt.ylim(np.min(pca[:, 0])-5, np.max(pca[:, 0])+5)
        plt.xlim(np.min(pca[:, 1])-5, np.max(pca[:, 1])+5)

    plt.tight_layout()
    # plt.show(block=True)
    plt.savefig(output_file_name + '.pdf', dpi=300, bbox_inches='tight')




def calcErrorRate(fastaFile, refSeq, bc_1, bc_2):

    correct_nt_array = np.zeros((len(refSeq)))
    total_nt_array = np.zeros((len(refSeq)))
    error_array = np.zeros((len(refSeq)))

    pat_start = r"(" + bc_1 + ")"
    pat_start = re.compile(pat_start)

    pat_end = r"(" + bc_2 + ")"
    pat_end = re.compile(pat_end)

    for header, seq in readFasta(fastaFile):

        search_start = pat_start.search(seq)
        search_end = pat_end.search(seq)

        if search_start and search_end:
            seq_trimmed = seq[search_start.span()[1]:search_end.span()[0]]
            if len(seq_trimmed) == len(refSeq):

                for i, nt in enumerate(seq_trimmed):
                    if nt == refSeq[i]:
                        correct_nt_array[i] += 1
                    total_nt_array[i] += 1


    for i in range(len(error_array)):
        error_array[i] = (1-(correct_nt_array[i]/total_nt_array[i]))*100.0

    print(correct_nt_array)
    print(total_nt_array)
    print(error_array)
    print("Mean error: %.3f%%" % (np.mean(error_array)))

    plt.xlabel("Position")
    plt.ylabel("Error rate (%)")
    plt.ylim(0, 100)
    plt.bar(np.arange(error_array.shape[0]), error_array)

    plt.tight_layout()
    # plt.show(block=True)
    plt.savefig('201215_randomer_fig3_error_rate_distribution.pdf', dpi=300, bbox_inches='tight')




#### Main #######

ref_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

# ref_seq = "GAAGGAAGAAAATGCAGAAAAAAAGAAAAAAATGTCTGG"

# T5 A1
file_list = ["T5_fasta/T5R0modNew.fasta", "T5_fasta/T5R1modNew.fasta", "T5_fasta/T5R2modNew.fasta",
"T5_fasta/T5R3modNew.fasta", "T5_fasta/T5R4modNew.fasta", "T5_fasta/T5R5modNew.fasta", "T5_fasta/T5R6modNew.fasta", "T5_fasta/T5R7modNew.fasta", "T5_fasta/T5R8_trim.fasta"]
bc_1 = "GGATGAGCGACGCTG"
bc_2 = "GTCTGGCGACTGCTC"
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=True, output_file_name="201229_T5_A1_HD_from_A1")
figTwoManifold(file_list, bc_1, bc_2, output_file_name="210115_randomer_fig2_tSNE_T5_A1_sequence_per_round_exact_motif",
    gen_random=True, labels=['R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', "Random"])

# ## T7 A2
file_list = ["demux/T7_R0.fasta", "demux/T7_R1.fasta", "demux/T7_R12.fasta", "demux/T7_R16.fasta"]
bc_1 = "GCACCGTGGACACAG"
bc_2 = "GACTGCCAGGTCGAG"
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=True, labels=["R0", "R1", "R12", "R16"], output_file_name="201229_T7_A2_HD_from_A1")
# figTwoManifold(file_list, bc_1, bc_2, output_file_name="201229_randomer_fig2_tSNE_T7_A2_sequence_per_round",
#     gen_random=True, labels=['R0', 'R1','R12', 'R16', "Random"])


# ## T8 A3
file_list = ["demux/T8_R0.fasta", "demux/T8_R1.fasta", "demux/T8_R12.fasta", "demux/T8_R16.fasta"]
bc_1 = "GCCTGCAAGTGCGAC"
bc_2 = "GAGACCACCGACGTG"
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=True, labels=["R0", "R1", "R12", "R16"], output_file_name="201229_T8_A3_HD_from_A1")
# figTwoManifold(file_list, bc_1, bc_2, output_file_name="201229_randomer_fig2_tSNE_T8_A3_sequence_per_round",
#     gen_random=True, labels=['R0', 'R1','R12', 'R16', "Random"])



# ref_seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

# ## T4 T1
file_list = ["demux/T4_R20.fasta", "demux/T4_R24.fasta"]
bc_1 = "GGCACTGCAGAGTCG"
bc_2 = "GACTCGAGCCAGAGC"
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=True, labels=["R20", "R24"], output_file_name="201229_T4_T1_HD_from_T1")
# figTwoManifold(file_list, bc_1, bc_2, output_file_name="201229_randomer_fig2_tSNE_T4_T1_sequence_per_round",
#     gen_random=True, labels=['R20', 'R24', "Random"])


# # ## T6 T2
file_list = ["demux/T6_R20.fasta", "demux/T6_R24.fasta"]
bc_1 = "GGACAGCGCGTGTAG"
bc_2 = "CTCCGTCGCTGCTAG"
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=True, labels=["R20", "R24"], output_file_name="201229_T6_T2_HD_from_T1")
# figTwoManifold(file_list, bc_1, bc_2, output_file_name="201229_randomer_fig2_tSNE_T6_T2_sequence_per_round",
#     gen_random=True, labels=['R20', 'R24', "Random"])



# plot_entropy_over_rounds(file_list, bc_1, bc_2)

# plot_hamming_distance_over_rounds(file_list, bc_1, bc_2, ref_seq)

# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, add_rand=True, write_csv=False)

# figTwoManifold(file_list, bc_1, bc_2)

# ref_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
# calcErrorRate("T5_fasta/T5R0modNew.fasta", ref_seq, bc_1, bc_2)




## GTP aptamer.

# file_list = ["demux/T8_R12.fasta", "demux/T8_R16.fasta"]
# bc_1 = "GCCTGCAAGTGCGAC"
# bc_2 = "GAGACC"

# ref_seq = "GAATGGGAAATAGCAAAAGAAAACAGAGACACCGACATGTTA"
# labels = ["R12", "R16"]
# plot_hamming_distance_over_rounds_stacked(file_list, bc_1, bc_2, ref_seq, labels)






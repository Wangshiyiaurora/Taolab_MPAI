import scrublet as scr
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

### Young_MCAI016
infile = "MCAI016_count.txt"
outfile = "MCAI016_Scrublet_result.txt"
finallist = []
with open(infile, 'r') as f:
    header = next(f)
    cell_barcodes = header.rstrip().split('\t')
    for line in f:
        tmpline = line.rstrip().split('\t')[1: ]
        tmplist = [float(s) for s in tmpline]
        finallist.append(tmplist)

finalarray = np.array(finallist)
count_matrix = np.transpose(finalarray)

scrub = scr.Scrublet(count_matrix, expected_doublet_rate = 0.08)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
scrub.plot_histogram()
plt.savefig('doublet_histogram.pdf', dpi = 300)
predicted_doublets_final = scrub.call_doublets(threshold = 0.25)

with open(outfile, 'w') as f:
    f.write('\t'.join(['CB', 'Scrublet', 'Scrublet_Score']) + '\n')
    for i in range(len(doublet_scores)):
        if predicted_doublets_final[i] == 0:
            result = 'Singlet'
        else:
            result = 'Doublet'
            f.write('\t'.join([cell_barcodes[i], result, str(doublet_scores[i])]) + '\n')

####tsne
def get_tsne(X, angle = 0.5, perplexity = 30, random_state = 0, verbose = False):
    from sklearn.manifold import TSNE
    return TSNE(angle = angle, perplexity = perplexity, random_state = random_state, verbose = verbose).fit_transform(X)
scrub.set_embedding('TSNE', scr.get_tsne(scrub.manifold_obs_, 0.5, 10))
scrub.plot_embedding('TSNE',order_points=True)
plt.savefig('Young_MCAI016/doublet_tsne.pdf', dpi = 300)

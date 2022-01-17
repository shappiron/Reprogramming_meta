import numpy as np
import multiprocessing
import string
import random
import pandas as pd
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

def GSEA(profile, geneset, return_info=False):
    if len(geneset) == 0:
        return np.nan
        
    ids = np.where(np.isin(profile, geneset))[0]
    ids_sort = np.sort(ids)

    n = geneset.shape[0]
    N = profile.shape[0]

    ids_sort += 1
    ids_stack = np.hstack([ids_sort, ids_sort-1, ids_sort+1])

    tile = np.tile(ids_sort, (3 * ids_sort.shape[0], 1))

    comp = np.sum(tile <= ids_stack[:, None], 1)
    phit =  comp / n
    pmiss = (ids_stack - comp) / (N - n)

    es = phit - pmiss
    id_max = np.argmax(np.abs(es))
    similarity_score = es[id_max]

    if return_info == True:
        ids_sort -= 1
        gene_names = profile[ids_sort]
        
        if similarity_score > 0:
            direction = "pos"
            significant = np.where(ids_sort <= id_max, "Y", "N")
            percent_top = np.sum(significant=='Y') / n * 100
            #percent_top <- sum(gene_table$Significant=="Y")/n*100
        else:
            direction = "neg"
            significant = np.where(ids_sort >= id_max, "Y", "N")
            percent_top = np.sum(significant=='Y') / n * 100

        info = {'direction':direction, 
                'top':percent_top,
                'genes':gene_names,
                'ids':ids_sort,
                'significant':significant,
                'es':es}
        return similarity_score, info
    else:
        return similarity_score

def premutation_biGSEA(profile_genes, termup, termdw):
    sample_up = np.random.choice(profile_genes, size=termup.shape[0], replace=False)
    sample_dw = np.random.choice(profile_genes, size=termdw.shape[0], replace=False)
    score_up = GSEA(profile_genes, sample_up, False)
    score_dw = GSEA(profile_genes, sample_dw, False)
    return score_up, score_dw


def association_test(profile, genesets_list, N_permutations=100, n_jobs=24, parallel=True):
    profile_name = profile.columns[0]
    profile_genes = profile.index.to_numpy()
    profile_sorted = profile.sort_values(profile_name, ascending=False)
    profile_sorted_genes = profile_sorted.index.to_numpy()

    out = pd.DataFrame(index=genesets_list.keys(), columns = ['Term', 'pval', 'padj', 'NES', 
                                                            'Top_up', 'Top_down', 
                                                            'Direction_up', 'Direction_down'])
    pool = multiprocessing.Pool(n_jobs)
    for term, double_set in tqdm(genesets_list.items()):
        termup = np.intersect1d(double_set['Up'], profile_genes)
        termdw = np.intersect1d(double_set['Down'], profile_genes)

        if (len(termup) + len(termdw)) == 0:
            continue
        
        if parallel:
            resampling_scores_up, resampling_scores_dw = zip(*pool.starmap(
                                                        premutation_biGSEA, 
                                                        [(profile_genes, termup, termdw)] * N_permutations))
        else:
            resampling_scores_up, resampling_scores_dw = zip(*[premutation_biGSEA(profile_genes, termup, termdw) 
                                                               for _ in range(N_permutations)])
        resampling_scores_up = np.asarray(resampling_scores_up)
        resampling_scores_dw = np.asarray(resampling_scores_dw)

        #normalize resampling scores by std
        sd_up = resampling_scores_up.std(ddof=1)
        sd_dw = resampling_scores_dw.std(ddof=1)
        resampling_scores_up = resampling_scores_up / sd_up
        resampling_scores_dw = resampling_scores_dw / sd_dw
        resampling_scores = (resampling_scores_up - resampling_scores_dw) / 2

        #run GSEA for profile
        termup = np.intersect1d(double_set['Up'], profile_sorted_genes)
        termdw = np.intersect1d(double_set['Down'], profile_sorted_genes)

        score_up, info_up = GSEA(profile_sorted_genes, termup, return_info=True)
        score_dw, info_dw = GSEA(profile_sorted_genes, termdw, return_info=True)
        info_up['es'] = info_up['es'] / sd_up
        info_dw['es'] = info_dw['es'] / sd_dw
        NES_up = score_up / sd_up
        NES_dw = score_dw / sd_dw

        final_score = (NES_up - NES_dw) / 2
        
        out.loc[term, "Term"] = term
        out.loc[term, "NES"] = final_score
        out.loc[term, "Direction_up"] = info_up['direction']
        out.loc[term, "Direction_down"] = info_dw['direction']
        out.loc[term, "Top_up"] = info_up['top']
        out.loc[term, "Top_down"] = info_dw['top']
        
    out['pval'] = out['NES'].apply(lambda x: np.sum(np.abs(resampling_scores) > np.abs(x)) / N_permutations)
    out['pval'] = np.where(out['pval']==0., 1/N_permutations, out['pval'])
    out['padj'] = multipletests(out['pval'], method='fdr_bh')[1]
    return out


def open_gmt_cmap(path_up, path_dw):
    with open(path_up) as f:
        gmtu = f.readlines()
    with open(path_dw) as f:    
        gmtd = f.readlines()
    terms = {}
    for s in gmtu:
        title, _, *genes  = s.strip().split('\t')
        terms[title] = {'Up':genes, 'Down':None}
    for s in gmtd:
        title, _, *genes  = s.strip().split('\t')
        terms[title]['Down'] = genes
    return terms  

#real data
if __name__ == '__main__':
    reprog = pd.read_csv('reprogramming_mouse_full.csv', index_col=0)
    profile = reprog[['symbol', 'logFC']]
    profile = profile.dropna(0)
    profile['symbol'] = profile['symbol'].apply(str.upper)
    profile = profile.set_index('symbol')

    path_up = 'drugs/cmap_up.gmt'
    path_dw = 'drugs/cmap_dw.gmt'
    genesets_list = open_gmt_cmap(path_up, path_dw)
    #
    #sub = {k:genesets_list[k] for k in list(genesets_list.keys())[0::200]}

    out = association_test(profile, genesets_list, 5000, 40, parallel=True)
    out.to_csv('bigsea_results_mouse_5000.csv')
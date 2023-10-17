import pandas as pd
import os
import scanpy as sc
import muon as mu
import numpy as np

# https://github.com/scverse/scanpy/issues/987
def obs_key_wise_subsampling(adata, obs_key, N):
    '''
    Subsample each class to same cell numbers (N). Classes are given by obs_key pointing to categorical in adata.obs.
    '''
    counts = adata.obs[obs_key].value_counts()
    # subsample indices per group defined by obs_key
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=True) for group in counts.index]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

figures = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/figures/rna'
sc_file = '/home/carsten/alvira_bioinformatics/postnatal_lung_multiome/data/single_cell_files'
deg_dir = f"{figures}/deg"
os.makedirs(deg_dir, exist_ok=True)
sc.set_figure_params(dpi=200, format="png")
sc.settings.figdir = figures

if __name__ == "__main__":
    mudata = mu.read(f'{sc_file}/multi_all_cells_processed.h5mu')
    adata = mudata.mod['rna']
    adata_norm = adata[adata.obs['treatment']=='Normoxia']
    sc.tl.rank_genes_groups(
        adata_norm,
        "lineage",
        method="wilcoxon",
        pts=True,
        key_added="rank_genes_groups_lineage",
    )
    with pd.ExcelWriter(
            f"{deg_dir}/lineage_markers.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs["lineage"].cat.categories:
            df = sc.get.rank_genes_groups_df(
                adata_norm, key="rank_genes_groups_lineage", group=ct
            )
            df.index = df['names']
            df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])


    print(adata)
    hyp_deg_dict = {}
    for lineage in adata.obs["lineage"].unique():
        print(lineage)
        figures_lin = f"{deg_dir}/{lineage}"
        os.makedirs(figures_lin, exist_ok=True)
        sc.settings.figdir = figures_lin
        lin_adata = adata[adata.obs["lineage"] == lineage]
        figures_lin_comp = f"{figures_lin}/cell_type_comparisons"
        os.makedirs(figures_lin_comp, exist_ok=True)
        ## broad ct markers
        for treat in lin_adata.obs['treatment'].cat.categories:
            treat_lin_adata = lin_adata[lin_adata.obs['treatment']==treat]
            subsampled = obs_key_wise_subsampling(treat_lin_adata, 'celltype',500)
            sc.tl.rank_genes_groups(
                subsampled,
                "celltype",
                method="wilcoxon",
                pts=True,
                key_added="rank_genes_groups_celltype",
            )
            with pd.ExcelWriter(
                    f"{figures_lin}/{lineage}_{treat.lower()}_celltype_markers.xlsx", engine="xlsxwriter"
            ) as writer:
                for ct in subsampled.obs["celltype"].cat.categories:
                    df= sc.get.rank_genes_groups_df(
                        subsampled, key="rank_genes_groups_celltype", group=ct
                    )
                    df.index = df['names']
                    df.to_excel(writer, sheet_name=f"{ct} v rest"[:31])

        ##  ct v every other ct
            figures_lin_treat_comp = f"{figures_lin}/cell_type_comparisons/{treat.lower()}"
            os.makedirs(figures_lin_treat_comp, exist_ok=True)
            for ct in treat_lin_adata.obs["celltype"].cat.categories:
                print(ct)
                with pd.ExcelWriter(
                        f"{figures_lin_treat_comp}/{ct}.xlsx", engine="xlsxwriter"
                ) as writer:
                    for ct2 in treat_lin_adata.obs["celltype"].cat.categories:
                        cts_adata = treat_lin_adata[treat_lin_adata.obs["celltype"].isin([ct, ct2])]
                        try:
                            sc.tl.rank_genes_groups(
                                cts_adata,
                                "celltype",
                                groups=[ct, ct2],
                                method="wilcoxon",
                                pts=True,
                                key_added="rank_genes_groups_celltype",
                            )
                            df = sc.get.rank_genes_groups_df(
                                cts_adata, key="rank_genes_groups_celltype", group=ct
                            )
                            df.index = df['names']
                            df.to_excel(writer, sheet_name=f"{ct} v {ct2}"[:31])
                        except:
                            print(treat)
                            print(ct)
                            print(ct2)
                            print('no comp')

        ## treatment
        with pd.ExcelWriter(
                f"{figures_lin}/hyperoxia_deg.xlsx", engine="xlsxwriter"
        ) as writer:
            for ct in lin_adata.obs["celltype"].cat.categories:
                ct_adata = lin_adata[lin_adata.obs["celltype"] == ct]
                try:
                    sc.tl.rank_genes_groups(
                        ct_adata,
                        "treatment",
                        method="wilcoxon",
                        pts=True,
                        key_added="rank_genes_groups_treatment",
                    )
                    df = sc.get.rank_genes_groups_df(
                        ct_adata,
                        key="rank_genes_groups_treatment",
                        group="Hyperoxia",
                    )
                    df.index = df['names']
                    df.to_excel(writer, sheet_name=f"{ct}")
                    hyp_deg_dict[ct]=df
                except:
                    print(ct)
                    print('no hyperoxia comparison')
                    continue
    gene_dict = {}
    for direction in ['up', 'down']:
        gene_dict[direction] = {}
        gene_dict[direction]['celltypes'] = {}
        for ct in adata.obs['celltype'].cat.categories:
            df = hyp_deg_dict[ct]
            df = df.loc[(df['pvals_adj'] < 0.05) |
                        (abs(df['scores']) > 1)]

            if direction == 'up':
                df = df.loc[df['logfoldchanges'] > 0]
            else:
                df = df.loc[df['logfoldchanges'] < 0]
            gene_ls = df.index.tolist()
            for gene in gene_ls:
                if gene in gene_dict[direction]['celltypes'].keys():
                    gene_dict[direction]['celltypes'][gene].append(ct)
                else:
                    gene_dict[direction]['celltypes'][gene] = [ct]

        gene_dict[direction]['number_cts'] = {}
        for gene in gene_dict[direction]['celltypes'].keys():
            gene_dict[direction]['number_cts'][gene] = len(gene_dict[direction]['celltypes'][gene])
    gene_dict_lin = {}
    for lineage in adata.obs['lineage'].cat.categories:
        gene_dict_lin[lineage] = {}
        lin_adata = adata[adata.obs['lineage'] == lineage]
        for direction in ['up', 'down']:
            gene_dict_lin[lineage][direction] = {}
            gene_dict_lin[lineage][direction]['celltypes'] = {}
            for ct in lin_adata.obs['celltype'].cat.categories:
                df = hyp_deg_dict[ct]
                df = df.loc[(df['pvals_adj'] < 0.05) |
                            (abs(df['scores']) > 1)]
                if direction == 'up':
                    df = df.loc[df['logfoldchanges'] > 0]
                else:
                    df = df.loc[df['logfoldchanges'] < 0]
                gene_ls = df.index.tolist()
                for gene in gene_ls:
                    if gene in gene_dict_lin[lineage][direction]['celltypes'].keys():
                        gene_dict_lin[lineage][direction]['celltypes'][gene].append(ct)
                    else:
                        gene_dict_lin[lineage][direction]['celltypes'][gene] = [ct]
            gene_dict_lin[lineage][direction]['number_cts'] = {}
            for gene in gene_dict_lin[lineage][direction]['celltypes'].keys():
                gene_dict_lin[lineage][direction]['number_cts'][gene] = len(
                    gene_dict_lin[lineage][direction]['celltypes'][gene])
    lin_ct_dict = {}
    for ct in adata.obs['celltype'].cat.categories:
        lin_ct_dict[ct] = adata[adata.obs['celltype'] == ct].obs['lineage'].unique()[0]
    with pd.ExcelWriter(
            f"{deg_dir}/hyperoxia_deg_curated.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs['celltype'].cat.categories:
            lin = lin_ct_dict[ct]
            df = hyp_deg_dict[ct]
            df = df.loc[(df['pvals_adj'] < 0.05) |
                        (abs(df['scores']) > 1)]
            df['number_ct_up'] = pd.Series({k: gene_dict['up']['number_cts'].get(k, 0) for k in df.index})
            df['number_ct_down'] = pd.Series({k: gene_dict['down']['number_cts'].get(k, 0) for k in df.index})
            df[f'number_{lin}_ct_up'] = pd.Series(
                {k: gene_dict_lin[lin]['up']['number_cts'].get(k, 0) for k in df.index})
            df[f'number_{lin}_ct_down'] = pd.Series(
                {k: gene_dict_lin[lin]['down']['number_cts'].get(k, 0) for k in df.index})
            df['cts_up'] = pd.Series({k: gene_dict['up']['celltypes'].get(k, 0) for k in df.index})
            df['cts_down'] = pd.Series({k: gene_dict['down']['celltypes'].get(k, 0) for k in df.index})
            df = df.sort_values('scores', ascending=False)
            df.to_excel(writer, sheet_name=ct[:31])
    with pd.ExcelWriter(
            f"{deg_dir}/hyperoxia_deg_all.xlsx", engine="xlsxwriter"
    ) as writer:
        for ct in adata.obs['celltype'].cat.categories:
            try:
                hyp_deg_dict[ct].to_excel(writer, sheet_name=f"{ct}")
            except:
                print(ct)
    print('DONE!')

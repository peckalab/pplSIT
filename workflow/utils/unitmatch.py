# # before loading data, update cluster_group.tsv
# # we have to "merge" the labels in cluster_KSlabel.tsv with the ones in cluster_group.tsv
# # cluster_group.tsv has priority, so only units that are not labeled in cluster_group.tsv will get labels from cluster_KSlabel.tsv
# for ks_dir in KS_dirs:
#     util.update_cluster_group_with_KSlabel(ks_dir)

def update_cluster_group_with_KSlabel(ks_dir):
    import os
    import pandas as pd

    cluster_group_path = os.path.join(ks_dir, 'cluster_group.tsv')
    cluster_KSlabel_path = os.path.join(ks_dir, 'cluster_KSLabel.tsv')

    if not os.path.isfile(cluster_group_path) or not os.path.isfile(cluster_KSlabel_path):
        print(f"Either {cluster_group_path} or {cluster_KSlabel_path} does not exist. Skipping update.")
        return

    # Read the existing cluster group and KS label files
    cluster_group_df = pd.read_csv(cluster_group_path, sep='\t')
    cluster_KSlabel_df = pd.read_csv(cluster_KSlabel_path, sep='\t')

    # if the headers of cluster_group.tsv are cluster_id	KSLabel then it is all done
    # if the headers are cluster_id	group then we need to update
    if 'group' in cluster_group_df.columns:
        # iterate over cluster_KSlabel_df, if the cluster_id is not in cluster_group_df, add it
        for _, row in cluster_KSlabel_df.iterrows():
            cluster_id = row['cluster_id']
            KSLabel = row['KSLabel']

            if cluster_id not in cluster_group_df['cluster_id'].values:
                # add new row
                new_row = pd.DataFrame({'cluster_id': [cluster_id], 'group': [KSLabel]})
                cluster_group_df = pd.concat([cluster_group_df, new_row], ignore_index=True)
    
        # now I just need to save the updated cluster_group_df back to cluster_group.tsv
        
        return cluster_group_df
    else:
        print(f"{cluster_group_path} already has KSLabel column. No update needed.")
        return None
        
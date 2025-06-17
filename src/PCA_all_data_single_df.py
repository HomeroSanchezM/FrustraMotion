import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from umap import UMAP
import argparse
import os
from sklearn.neighbors import NearestNeighbors
from community import community_louvain  # python-louvain package
import networkx as nx

monomer_id = {
    '0':1, '1':2, '2':3, '3':4, '4':5, '5':6, '6':7, '7':8, 'A':9, 'B':10,
    'C':11,'D':12,'E':13,'F':14,'G':15,'H':16,'I':17,'J':18,'K':19,'L':20,
    'M':21,'N':22,'O':23,'P':24,'Q':25,'R':26,'S':27,'T':28,'U':29,'V':30,
    'W':31,'X':32,'Y':33,'Z':34,'a':35,'b':36,'c':37,'d':38,'e':39,'f':40,
    'g':41,'h':42,'i':43,'j':44,'k':45,'l':46,'m':47,'n':48,'o':49,'p':50,
    'q':51,'r':52,'s':53,'t':54,'u':55,'v':56,'w':57,'x':58,'y':59,'z':60
}

# inverse dico
id_monomer = {v: k for k, v in monomer_id.items()}

def parse_arguments():
    """Parse command line arguments with flexible options"""
    parser = argparse.ArgumentParser(description='Advanced PCA and UMAP analysis')
    parser.add_argument('matrix_file', help='Path to the precomputed matrix file')
    parser.add_argument('--pca_components', type=int, default=60,
                        help='Number of PCA components to compute')
    parser.add_argument('--umap_components', type=int, choices=[6, 10, 20, 50, 60],
                        default=20, help='Top PCs to use for UMAP')
    parser.add_argument('--umap_neighbors', type=int, default=15,
                        help='UMAP n_neighbors parameter')
    parser.add_argument('--umap_min_dist', type=float, default=0.1,
                        help='UMAP min_dist parameter')
    parser.add_argument('--louvain_resolution', type=float, default=1.0,
                        help='Resolution parameter for Louvain clustering')
    return parser.parse_args()


def plot_elbow(pca, output_prefix):
    """Enhanced elbow plot with additional information"""
    plt.figure(figsize=(18, 6))
    plt.plot(range(1, len(pca.explained_variance_ratio_) + 1),
             pca.explained_variance_ratio_, 'o-', linewidth=2, label='Individual')
    plt.xlabel('Number of Components')
    plt.ylabel('Explained Variance Ratio')
    plt.title(f'{output_prefix[:2]} Elbow Plot (Total Variance: {np.sum(pca.explained_variance_ratio_):.2f})')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xticks(range(0, len(pca.explained_variance_ratio_) + 1, 5))
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_elbow_plot.png', dpi=300)
    plt.show()


def plot_projection(embedding, title, output_file, monomer_groups=None, pca_var=None, clusters=None):
    """Universal plotting function for both PCA and UMAP with optional clustering"""
    plt.figure(figsize=(12, 12))

    if clusters is not None:
        # Plot with Louvain clusters
        unique_clusters = np.unique(clusters)
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))

        for cluster, color in zip(unique_clusters, colors):
            mask = clusters == cluster
            plt.scatter(embedding[mask, 0], embedding[mask, 1],
                        s=150, color=color, label=f'Cluster {cluster}')

            # Add labels for points in this cluster
            for i in np.where(mask)[0]:
                plt.text(embedding[i, 0], embedding[i, 1], f'M{i + 1}',
                         fontsize=9, ha='center', va='bottom')

        plt.legend(title="Louvain Clusters", bbox_to_anchor=(1.05, 1))

    elif monomer_groups is not None:
        # Original plotting with monomer groups
        colors = plt.cm.tab20(np.linspace(0, 1, len(monomer_groups)))
        group_colors = {group: colors[i] for i, group in enumerate(monomer_groups.keys())}

        monomer_to_group = {}
        for group, monomers in monomer_groups.items():
            for m in monomers:
                monomer_to_group[m] = group

        for i in range(embedding.shape[0]):
            group = monomer_to_group.get(i + 1, 'other')
            color = group_colors.get(group, 'gray')
            plt.scatter(embedding[i, 0], embedding[i, 1],
                        s=150, color=color, label=f'M{i + 1} ({group})')
            plt.text(embedding[i, 0], embedding[i, 1], f'M{i + 1}',
                     fontsize=9, ha='center', va='bottom')

        handles = [plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=group_colors[g], markersize=10)
                   for g in monomer_groups]

    # Add axis labels
    if pca_var is not None:
        plt.xlabel(f'PC1 ({pca_var[0] * 100:.1f}%)')
        plt.ylabel(f'PC2 ({pca_var[1] * 100:.1f}%)')
    else:
        plt.xlabel(f'{title.split()[0]} Dimension 1')
        plt.ylabel(f'{title.split()[0]} Dimension 2')

    plt.title(title)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()


def perform_louvain_clustering(embedding, n_neighbors=15, resolution=1.0):
    """Perform Louvain clustering on UMAP embedding"""
    # Create k-nearest neighbors graph
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(embedding)
    distances, indices = nbrs.kneighbors(embedding)

    # Create a weighted graph
    G = nx.Graph()

    for i in range(len(embedding)):
        for j, dist in zip(indices[i], distances[i]):
            if i != j:  # no self-loops
                # Convert distance to similarity (weight)
                weight = 1.0 / (1.0 + dist)
                G.add_edge(i, j, weight=weight)

    # Perform Louvain clustering
    partition = community_louvain.best_partition(G, resolution=resolution, random_state=42)

    # Convert partition to cluster labels
    clusters = np.array([partition[i] for i in range(len(embedding))])

    return clusters


def main():
    global args
    args = parse_arguments()

    # 1. Load precomputed matrix
    print(f"\nLoading precomputed matrix from {args.matrix_file}")
    df = pd.read_csv(args.matrix_file, sep='\t', header=None)

    # Clean the data (remove NaN columns and rows)
    df = df.dropna(axis=1, how='all')
    df = df.dropna(axis=0, how='any')

    print(f"\nMatrix shape: {df.shape}")
    print("Sample data:")
    print(df.iloc[:5, :5].to_string())

    # 2. Standardization
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df)

    # 3. PCA Analysis
    output_prefix = os.path.splitext(os.path.basename(args.matrix_file))[0]

    pca = PCA(n_components=args.pca_components)
    pca_components = pca.fit_transform(scaled_data)

    # Plot elbow
    plot_elbow(pca, output_prefix)

    # 4. UMAP Analysis
    print(f"\nRunning UMAP with parameters:")
    print(f"- Using top {args.umap_components} PCs")
    print(f"- n_neighbors: {args.umap_neighbors}")
    print(f"- min_dist: {args.umap_min_dist}")

    reducer = UMAP(n_components=2,
                   n_neighbors=args.umap_neighbors,
                   min_dist=args.umap_min_dist,
                   random_state=42)

    umap_embedding = reducer.fit_transform(pca_components[:, :args.umap_components])

    # 5. Louvain Clustering
    
    print(f"\nPerforming Louvain clustering with resolution={args.louvain_resolution}")
    clusters = perform_louvain_clustering(umap_embedding,
                                          n_neighbors=args.umap_neighbors,
                                          resolution=args.louvain_resolution)

    print("\nCluster assignments with monomer names:")
    cluster_dict = {}
    for i, cluster in enumerate(clusters):
        monomer_num = i + 1  # Les monomères sont numérotés de 1 à 60
        monomer_name = id_monomer.get(monomer_num, f"Unknown({monomer_num})")
        print(f"M{monomer_num} ({monomer_name}): Cluster {cluster}")

        # Stocker les résultats par cluster
        if cluster not in cluster_dict:
            cluster_dict[cluster] = []
        cluster_dict[cluster].append((monomer_num, monomer_name))

    # Afficher la liste regroupée par cluster
    print("\nMonomer groups by cluster:")
    for cluster in sorted(cluster_dict.keys()):
        monomers = cluster_dict[cluster]
        print(f"\nCluster {cluster}:")
        print("Numbers:", ", ".join(f"M{m[0]}" for m in monomers))
        print("Names:  ", ", ".join(m[1] for m in monomers))
    print("\nCluster assignments:")
    for i, cluster in enumerate(clusters):
        print(f"M{i + 1}: Cluster {cluster}")

    # 6. Visualization
    monomer_groups = {
        'group1': (1, 6, 46, 51, 56),
        'group2': (14, 15, 16, 17, 18),
        'group3': (2, 7, 47, 52, 57),
        'group4': (3, 8, 48, 53, 58),
        'group5': (4, 44, 49, 54, 59),
        'group6': (5, 45, 50, 55, 60),
        'group7': (9, 10, 11, 12, 13),
        'group8': (19, 24, 29, 34, 39),
        'group9': (20, 25, 30, 35, 40),
        'group10': (21, 26, 31, 36, 41),
        'group11': (22, 27, 32, 37, 42),
        'group12': (23, 28, 33, 38, 43)
    }

    # Plot PCA (first 2 components) with monomer groups
    plot_projection(pca_components[:, :2],
                    f' {output_prefix[:2]} PCA of Frustration Landscape (First 2 PCs)\nVariance: {np.sum(pca.explained_variance_ratio_[:2]):.2f}',
                    f"{output_prefix}_pca_plot.png",
                    monomer_groups,
                    pca_var=pca.explained_variance_ratio_[:2])

    # Plot UMAP with monomer groups
    plot_projection(umap_embedding,
                    f' {output_prefix[:2]} UMAP of Frustration Landscape (Top {args.umap_components} PCs)\nParams: n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}',
                    f"{output_prefix}_umap_plot_monomer_groups.png",
                    monomer_groups)

    # Plot UMAP with Louvain clusters
    plot_projection(umap_embedding,
                    f' {output_prefix[:2]} UMAP with Louvain Clustering (Resolution={args.louvain_resolution})\nTop {args.umap_components} PCs, n_neighbors={args.umap_neighbors}',
                    f"{output_prefix}_umap_plot_louvain_clusters.png",
                    clusters=clusters)


if __name__ == "__main__":
    main()
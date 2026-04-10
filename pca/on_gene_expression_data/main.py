import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

def get_gata3_id(filepath):
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            cols = line.split('\t')
            if len(cols) > 4 and cols[4].strip() == 'GATA3':
                return cols[0].strip()
    return "4359"

def main():
    labels_df = pd.read_csv('data/class.tsv', sep='\t', header=None)
    y = labels_df.iloc[:, 0].values

    xbp1_id = "4404"
    gata3_id = get_gata3_id('data/columns.tsv.gz')

    df = pd.read_csv('data/filtered.tsv.gz', sep='\t')
    df.columns = df.columns.astype(str).str.strip()

    X = df[[gata3_id, xbp1_id]].copy()
    X.columns = ['GATA3', 'XBP1']
    X['Class'] = y
    X['Label'] = X['Class'].map({1: 'ER+', 0: 'ER-'})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    palette = {'ER+': '#d62728', 'ER-': '#1f77b4'}

    sns.scatterplot(data=X, x='GATA3', y='XBP1', hue='Label', palette=palette, ax=ax1)
    ax1.set_title('Gene Expression')
    ax1.set_xlabel('GATA3 Expression')
    ax1.set_ylabel('XBP1 Expression')

    pca = PCA(n_components=1)
    X['PC1'] = pca.fit_transform(X[['GATA3', 'XBP1']])

    sns.scatterplot(data=X, x='PC1', y=2, hue='Label', palette=palette, ax=ax2)
    sns.scatterplot(data=X[X['Class'] == 1], x='PC1', y=1, color='#d62728', ax=ax2)
    sns.scatterplot(data=X[X['Class'] == 0], x='PC1', y=0, color='#1f77b4', ax=ax2)

    ax2.set_title('PC1 Projection')
    ax2.set_xlabel('PC1')
    ax2.set_yticks([0, 1, 2])
    ax2.set_yticklabels(['ER-', 'ER+', 'All'])

    plt.tight_layout()
    plt.savefig('pca_repro.png')

if __name__ == "__main__":
    main()

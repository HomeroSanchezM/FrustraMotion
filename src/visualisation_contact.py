import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from collections import defaultdict
import numpy as np

# Configuration des frames par défaut à visualiser
DEFAULT_FRAMES = [0, 4900]  # Modifier cette liste pour changer les frames par défaut


def load_contact_data(protein_name, isolate_mode):
    """Charge les données de contacts précédemment sauvegardées"""
    data_dir = os.path.join('../contact_data', 'Isolated' if isolate_mode else 'Not_isolated', protein_name)

    if not os.path.exists(data_dir):
        print(f"Error: No data found for protein {protein_name}")
        return None

    contact_files = [f for f in os.listdir(data_dir) if f.endswith('.tsv')]
    if not contact_files:
        print(f"Error: No TSV files found in {data_dir}")
        return None

    if len(contact_files) > 1:
        # Multiples chaînes - chargement dans un dictionnaire
        data = {}
        for f in contact_files:
            chain_id = f.split('_')[-1].split('.')[0]
            df = pd.read_csv(os.path.join(data_dir, f), sep='\t')
            # Conversion des Frames en numérique et tri
            df['Frame'] = pd.to_numeric(df['Frame'], errors='coerce')
            df = df.sort_values('Frame')
            data[chain_id] = df
        return data
    else:
        # Fichier unique
        df = pd.read_csv(os.path.join(data_dir, contact_files[0]), sep='\t')
        df['Frame'] = pd.to_numeric(df['Frame'], errors='coerce')
        return df.sort_values('Frame')


def filter_by_frames(df, frames):
    """Filtre le dataframe pour ne garder que les frames spécifiées"""
    return df[df['Frame'].isin(frames)]


def plot_contact_frustration(contact_data, output_prefix):
    """Génère différents graphiques à partir des données de contacts"""
    if isinstance(contact_data, dict):
        # Multiples chaînes - graphique pour chaque chaîne
        for chain_id, df in contact_data.items():
            _plot_single_chain(df, f"{output_prefix}_chain_{chain_id}")
    else:
        # Dataframe unique
        _plot_single_chain(contact_data, output_prefix)


def _plot_single_chain(df, output_prefix):
    """Génère les graphiques pour une seule chaîne"""
    # 1. Distribution des frustrations
    plt.figure(figsize=(10, 6))
    sns.histplot(df['FrstIndex'], bins=30, kde=True)
    plt.title('Distribution des indices de frustration')
    plt.xlabel('Indice de frustration')
    plt.ylabel('Nombre')
    plt.savefig(f"{output_prefix}_frustration_dist.png")
    plt.close()

    # 2. Frustration par type de puit
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Welltype', y='FrstIndex', data=df)
    plt.title('Frustration par type de puit')
    plt.savefig(f"{output_prefix}_frustration_by_welltype.png")
    plt.close()

    # 3. Proportions des états de frustration
    plt.figure(figsize=(8, 8))
    df['FrstState'].value_counts().plot.pie(autopct='%1.1f%%')
    plt.title('Proportion des états de frustration')
    plt.savefig(f"{output_prefix}_frustration_states.png")
    plt.close()

    # 4. Évolution de la frustration moyenne sur les frames
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        x='Frame', y='FrstIndex',
        data=df.groupby('Frame')['FrstIndex'].mean().reset_index(),
        marker='o'
    )
    plt.title('Évolution moyenne de la frustration')
    plt.xlabel('Numéro de frame')
    plt.ylabel('Indice moyen de frustration')
    plt.grid(True)
    plt.savefig(f"{output_prefix}_frustration_evolution.png")
    plt.close()


def analyze_residue_contacts(contact_data, residue_id):
    """Analyse tous les contacts pour un résidu spécifique"""
    if isinstance(contact_data, dict):
        # Trouve la chaîne du résidu
        chain_id = residue_id.split(':')[0]
        if chain_id not in contact_data:
            print(f"Error: Chaîne {chain_id} non trouvée dans les données")
            return None, None
        df = contact_data[chain_id]
    else:
        df = contact_data

    # Tous les contacts impliquant ce résidu
    contacts = df[(df['ResID1'] == residue_id) | (df['ResID2'] == residue_id)].copy()

    # Ajoute l'information de direction
    contacts['Direction'] = contacts.apply(
        lambda row: 'sortant' if row['ResID1'] == residue_id else 'entrant',
        axis=1
    )

    # Calcule les statistiques
    stats = {
        'total_contacts': len(contacts),
        'avg_frustration': contacts['FrstIndex'].mean(),
        'highly_frustrated': sum(contacts['FrstState'] == 'highly'),
        'minimally_frustrated': sum(contacts['FrstState'] == 'minimally'),
        'neutral': sum(contacts['FrstState'] == 'neutral'),
        'frames_analyzed': contacts['Frame'].nunique(),
        'first_frame': contacts['Frame'].min(),
        'last_frame': contacts['Frame'].max()
    }

    return contacts, stats


def plot_residue_evolution(contacts, residue_id, output_prefix):
    """Trace l'évolution de la frustration pour un résidu spécifique"""
    if contacts.empty:
        return

    plt.figure(figsize=(14, 7))

    # Contacts sortants
    outgoing = contacts[contacts['Direction'] == 'sortant']
    if not outgoing.empty:
        sns.lineplot(
            x='Frame', y='FrstIndex',
            data=outgoing.groupby('Frame')['FrstIndex'].mean().reset_index(),
            label='Contacts sortants',
            color='red',
            marker='o'
        )

    # Contacts entrants
    incoming = contacts[contacts['Direction'] == 'entrant']
    if not incoming.empty:
        sns.lineplot(
            x='Frame', y='FrstIndex',
            data=incoming.groupby('Frame')['FrstIndex'].mean().reset_index(),
            label='Contacts entrants',
            color='blue',
            marker='o'
        )

    plt.title(f'Évolution de la frustration pour le résidu {residue_id}')
    plt.xlabel('Numéro de frame')
    plt.ylabel('Indice de frustration')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{output_prefix}_residue_evolution.png")
    plt.close()


def plot_residue_contacts(contacts, residue_id, output_prefix, frames=None):
    """Visualisation détaillée des contacts pour des frames spécifiques ou toutes les frames"""
    if contacts.empty:
        return None, None

    # Si aucun frame spécifique n'est donné, on prend toutes les frames disponibles
    if frames is None:
        frames = sorted(contacts['Frame'].unique())
        frame_label = "all_frames"
    else:
        contacts = filter_by_frames(contacts, frames)
        frame_label = f"frames_{'_'.join(map(str, frames))}"

    # Création d'une palette de couleurs pour les différents contacts
    unique_pairs = contacts['ResID1'] + "-" + contacts['ResID2']
    unique_pairs = unique_pairs.unique()
    palette = sns.color_palette("husl", len(unique_pairs))
    color_map = {pair: palette[i] for i, pair in enumerate(unique_pairs)}

    # Mapping des tailles pour les différents types de contacts
    size_map = {
        'short': 100,
        'long': 60,
        'water-mediated': 40
    }

    plt.figure(figsize=(20, 10))

    # Tracer chaque paire de contacts séparément
    for pair in unique_pairs:
        pair_contacts = contacts[(contacts['ResID1'] + "-" + contacts['ResID2'] == pair) |
                                 (contacts['ResID2'] + "-" + contacts['ResID1'] == pair)]

        # Trier par frame pour une bonne connexion des points
        pair_contacts = pair_contacts.sort_values('Frame')

        # Tracer les points
        sns.scatterplot(
            x='Frame',
            y='FrstIndex',
            data=pair_contacts,
            color=color_map[pair],
            s=pair_contacts['Welltype'].map(size_map),
            alpha=0.8,
            label=pair
        )

        # Connecter les points du même contact avec une ligne
        plt.plot(
            pair_contacts['Frame'],
            pair_contacts['FrstIndex'],
            color=color_map[pair],
            linestyle='-',
            alpha=0.4,
            linewidth=1
        )

    # Inversion de l'axe Y
    plt.gca().invert_yaxis()

    # Ajouter des lignes horizontales pointillées
    plt.axhline(y=-1, color='red', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(y=0.78, color='green', linestyle='--', linewidth=1, alpha=0.7)

    # Configuration du graphique
    plt.title(f'Évolution des contacts pour le résidu {residue_id}\nFrames: 0-4900')
    plt.xlabel('Numéro de frame')
    plt.ylabel('Indice de frustration (inversé)')
    plt.grid(True)

    # Amélioration de l'affichage des ticks sur l'axe X
    plt.xticks(np.arange(0, 4901, 500))

    # Légende personnalisée
    handles = []
    for pair in unique_pairs:
        handles.append(plt.Line2D([0], [0],
                                  marker='o',
                                  color='w',
                                  markerfacecolor=color_map[pair],
                                  markersize=10,
                                  label=pair))

    # Légende pour les types de contacts
    size_legend = [
        plt.Line2D([0], [0],
                   marker='o',
                   color='w',
                   markerfacecolor='gray',
                   markersize=np.sqrt(size_map['short']) / 2,
                   label='Short'),
        plt.Line2D([0], [0],
                   marker='o',
                   color='w',
                   markerfacecolor='gray',
                   markersize=np.sqrt(size_map['long']) / 2,
                   label='Long'),
        plt.Line2D([0], [0],
                   marker='o',
                   color='w',
                   markerfacecolor='gray',
                   markersize=np.sqrt(size_map['water-mediated']) / 2,
                   label='Water-mediated')
    ]

    # Première légende pour les paires de contacts
    first_legend = plt.legend(handles=handles,
                              title='Paires de contacts',
                              bbox_to_anchor=(1.05, 1),
                              loc='upper left')

    # Ajouter la légende des tailles
    plt.gca().add_artist(first_legend)
    plt.legend(handles=size_legend,
               title='Type de contact',
               bbox_to_anchor=(1.05, 0.7),
               loc='upper left')

    plt.tight_layout()

    # Obtenir les limites de l'axe Y avant de sauvegarder
    y_limits = plt.gca().get_ylim()

    plt.savefig(f"{output_prefix}_contact_evolution_{frame_label}.png",
                bbox_inches='tight',
                dpi=300)
    plt.close()

    return y_limits, color_map


def plot_individual_contacts(contacts, residue_id, output_prefix, color_map, global_y_limits):
    """Génère des graphiques individuels pour chaque contact"""
    if contacts.empty:
        return

    # Créer un sous-répertoire pour les graphiques individuels
    individual_dir = os.path.join(os.path.dirname(output_prefix), "individual_contacts")
    os.makedirs(individual_dir, exist_ok=True)

    # Mapping des tailles pour les différents types de contacts
    size_map = {
        'short': 100,
        'long': 60,
        'water-mediated': 40
    }

    # Obtenir toutes les paires de contacts uniques
    unique_pairs = contacts['ResID1'] + "-" + contacts['ResID2']
    unique_pairs = unique_pairs.unique()

    for pair in unique_pairs:
        pair_contacts = contacts[(contacts['ResID1'] + "-" + contacts['ResID2'] == pair) |
                                 (contacts['ResID2'] + "-" + contacts['ResID1'] == pair)]

        # Trier par frame
        pair_contacts = pair_contacts.sort_values('Frame')

        plt.figure(figsize=(12, 6))

        # Tracer les points
        sns.scatterplot(
            x='Frame',
            y='FrstIndex',
            data=pair_contacts,
            color=color_map[pair],
            s=pair_contacts['Welltype'].map(size_map),
            alpha=0.8
        )

        # Connecter les points avec une ligne
        plt.plot(
            pair_contacts['Frame'],
            pair_contacts['FrstIndex'],
            color=color_map[pair],
            linestyle='-',
            alpha=0.6,
            linewidth=1.5
        )

        # Utiliser les mêmes limites d'axe Y que le graphique global
        plt.ylim(global_y_limits)

        # Inversion de l'axe Y et X
        #plt.gca().invert_yaxis()
        #plt.gca().invert_xaxis()

        # Ajouter des lignes horizontales pointillées
        plt.axhline(y=-1, color='red', linestyle='--', linewidth=1, alpha=0.7)
        plt.axhline(y=0.78, color='green', linestyle='--', linewidth=1, alpha=0.7)

        # Configuration du graphique
        plt.title(f'Évolution du contact {pair}\npour le résidu {residue_id}')
        plt.xlabel('Numéro de frame (inversé)')
        plt.ylabel('Indice de frustration (inversé)')
        plt.grid(True)

        # Amélioration de l'affichage des ticks sur l'axe X
        plt.xticks(np.arange(0, 4901, 500))

        plt.tight_layout()

        # Sauvegarder le graphique
        pair_name = pair.replace(':', '_').replace('-', '_')
        plt.savefig(
            os.path.join(individual_dir, f"{os.path.basename(output_prefix)}_{pair_name}_individual.png"),
            bbox_inches='tight',
            dpi=300
        )
        plt.close()
        
def main():
    parser = argparse.ArgumentParser(description='Analyse des données de frustration des contacts.')
    parser.add_argument('protein', type=str, help='Nom de la protéine à analyser')
    parser.add_argument('--isolate', action='store_true', help='Analyser les chaînes isolées')
    parser.add_argument('--residue', type=str, help='Résidu spécifique à analyser (format: Chain:AAResNum)')
    parser.add_argument('--frames', nargs='+', type=int, help='Frames spécifiques à analyser (ex: 0 10)')
    parser.add_argument('--plot_all', action='store_true', help='Générer tous les graphiques disponibles')

    args = parser.parse_args()

    # Utilise les frames fournies ou celles par défaut
    frames_to_plot = args.frames if args.frames else DEFAULT_FRAMES

    # Charge les données
    contact_data = load_contact_data(args.protein, args.isolate)
    if contact_data is None:
        return

    # Crée le répertoire de sortie
    output_dir = os.path.join('../contact_analysis', args.protein)
    os.makedirs(output_dir, exist_ok=True)

    # Génère les graphiques de base
    plot_contact_frustration(contact_data, os.path.join(output_dir, args.protein))

    # Analyse un résidu spécifique si demandé
    if args.residue:
        contacts, stats = analyze_residue_contacts(contact_data, args.residue)

        if contacts is not None:
            print(f"\nAnalyse pour le résidu {args.residue}:")
            print(f"Contacts totaux: {stats['total_contacts']}")
            print(f"Frustration moyenne: {stats['avg_frustration']:.3f}")
            print(f"Contacts très frustrés: {stats['highly_frustrated']}")
            print(f"Contacts peu frustrés: {stats['minimally_frustrated']}")
            print(f"Contacts neutres: {stats['neutral']}")
            print(f"Frames analysées: {stats['frames_analyzed']} (de {stats['first_frame']} à {stats['last_frame']})")

            # Sauvegarde les contacts dans un fichier
            residue_file = os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}.tsv")
            contacts.to_csv(residue_file, sep='\t', index=False)
            print(f"\nContacts détaillés sauvegardés dans {residue_file}")

            # Trace l'évolution sur toutes les frames
            y_limits, color_map = plot_residue_contacts(
                contacts, args.residue,
                os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}"),
                frames=args.frames if args.frames else None  # None pour toutes les frames
            )
            print(f"Graphique d'évolution sauvegardé pour les frames: {args.frames if args.frames else 'tous'}")

            # Si on a tracé toutes les frames, génère aussi les graphiques individuels
            if args.frames is None and y_limits is not None:
                plot_individual_contacts(
                    contacts, args.residue,
                    os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}"),
                    color_map, y_limits
                )
                print("Graphiques individuels des contacts sauvegardés dans le dossier 'individual_contacts'")

    # Génère des graphiques supplémentaires si demandé
    if args.plot_all:
        print("\nGénération de graphiques supplémentaires...")
        _generate_additional_plots(contact_data, os.path.join(output_dir, args.protein))


def _generate_additional_plots(contact_data, output_prefix):
    """Génère des graphiques supplémentaires"""
    if isinstance(contact_data, dict):
        for chain_id, df in contact_data.items():
            _plot_contact_heatmap(df, f"{output_prefix}_chain_{chain_id}")
    else:
        _plot_contact_heatmap(contact_data, output_prefix)


def _plot_contact_heatmap(df, output_prefix):
    """Heatmap des contacts à travers les frames"""
    try:
        # Crée une table pivot pour la heatmap
        pivot_df = df.pivot_table(
            index='ResID1',
            columns='Frame',
            values='FrstIndex',
            aggfunc='mean'
        )

        plt.figure(figsize=(15, 10))
        sns.heatmap(
            pivot_df,
            cmap='coolwarm',
            center=0,
            cbar_kws={'label': 'Indice de frustration'}
        )
        plt.title('Heatmap des frustrations des contacts')
        plt.xlabel('Numéro de frame')
        plt.ylabel('Résidu')
        plt.savefig(f"{output_prefix}_heatmap.png")
        plt.close()
    except Exception as e:
        print(f"Erreur lors de la génération de la heatmap: {e}")


if __name__ == "__main__":
    main()
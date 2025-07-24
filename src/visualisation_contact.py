import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from collections import defaultdict
import numpy as np

# Configuration des frames par défaut à visualiser
DEFAULT_FRAMES = [0, 4900]  # Modifier cette liste pour changer les frames par défaut


def load_contact_data(dataframes_dir):
    """Charge les données de contacts à partir du répertoire spécifié"""
    if not os.path.exists(dataframes_dir):
        print(f"Error: Directory {dataframes_dir} does not exist")
        return None

    contact_files = [f for f in os.listdir(dataframes_dir) if f.endswith('.csv')]
    if not contact_files:
        print(f"Error: No CSV files found in {dataframes_dir}")
        return None

    # Charger tous les fichiers et les combiner
    dfs = []
    for f in contact_files:
        df = pd.read_csv(os.path.join(dataframes_dir, f))
        df['Frame'] = pd.to_numeric(df['Frame'], errors='coerce')
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df.sort_values('Frame')


def filter_contacts(contacts, residue_id, contact_type=None):
    """Filtre les contacts selon le type demandé (inter/intra)"""
    # Tous les contacts impliquant ce résidu
    filtered = contacts[(contacts['ResID1'] == residue_id) | (contacts['ResID2'] == residue_id)].copy()

    if contact_type == 'intra':
        # Seulement les contacts intra-chaîne
        filtered = filtered[filtered['ChainRes1'] == filtered['ChainRes2']]
    elif contact_type == 'inter':
        # Seulement les contacts inter-chaînes
        filtered = filtered[filtered['ChainRes1'] != filtered['ChainRes2']]

    return filtered


def analyze_residue_contacts(contacts, residue_id):
    """Analyse tous les contacts pour un résidu spécifique"""
    if contacts.empty:
        return None, None

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
        'last_frame': contacts['Frame'].max(),
        'intra_chain': sum(contacts['ChainRes1'] == contacts['ChainRes2']),
        'inter_chain': sum(contacts['ChainRes1'] != contacts['ChainRes2'])
    }

    return contacts, stats


def plot_residue_contacts(contacts, residue_id, output_prefix,dataframe_dir,contact_type , frames ):
    """Visualisation des contacts pour un résidu spécifique"""
    # extract Protein_name and type
    protein_name = dataframe_dir.split('/')[-2]
    Type = dataframe_dir.split('/')[-3]
    if contacts.empty:
        return None, None

    # Si aucun frame spécifique n'est donné, on prend toutes les frames disponibles
    if frames is None:
        frames = sorted(contacts['Frame'].unique())
        frame_label = "all_frames"
    else:
        contacts = contacts[contacts['Frame'].isin(frames)]
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
    plt.title(f'Contacts frustration for residue {residue_id}\n ( {protein_name}, {Type}, Frames: {frame_label}, only-{contact_type} ) ')
    plt.xlabel('Numéro de frame')
    plt.ylabel('Indice de frustration (inversé)')
    plt.grid(True)

    # Amélioration de l'affichage des ticks sur l'axe X
    plt.xticks(np.arange(min(frames), max(frames) + 1, 500))

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

    plt.tight_layout()

    # Obtenir les limites de l'axe Y avant de sauvegarder
    y_limits = plt.gca().get_ylim()

    plt.savefig(f"{output_prefix}_contact_evolution_{frame_label}.png",
                bbox_inches='tight',
                dpi=300)
    plt.show()
    plt.close()

    return y_limits, color_map


def plot_individual_contacts(contacts, residue_id, output_prefix, color_map, global_y_limits):
    """Génère des graphiques individuels pour chaque contact avec mise en évidence d'un frame spécifique"""
    HIGHLIGHT_FRAME = 6000

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
        ax = plt.gca()

        # Tracer les points
        sns.scatterplot(
            x='Frame',
            y='FrstIndex',
            data=pair_contacts,
            color=color_map[pair],
            s=pair_contacts['Welltype'].map(size_map),
            alpha=0.8,
            ax=ax
        )

        # Connecter les points avec une ligne
        ax.plot(
            pair_contacts['Frame'],
            pair_contacts['FrstIndex'],
            color=color_map[pair],
            linestyle='-',
            alpha=0.6,
            linewidth=1.5
        )

        # Mettre en évidence le frame spécifié s'il existe
        if HIGHLIGHT_FRAME in pair_contacts['Frame'].values:
            highlight_data = pair_contacts[pair_contacts['Frame'] == HIGHLIGHT_FRAME].iloc[0]
            frustration = highlight_data['FrstIndex']

            # Déterminer la couleur en fonction de la frustration
            if frustration < -1:
                highlight_color = 'red'  # Très frustré
            elif frustration > 0.78:
                highlight_color = 'green'  # Peu frustré
            else:
                highlight_color = 'gray'  # Neutre

            # Mettre en évidence le point
            ax.scatter(
                x=HIGHLIGHT_FRAME,
                y=frustration,
                color=highlight_color,
                s=size_map[highlight_data['Welltype']] * 1.5,  # 50% plus grand
                alpha=1.0,
                zorder=4,
                edgecolor='black',
                linewidth=1
            )

            # Ajouter une zone verticale de mise en évidence
            ax.axvspan(
                HIGHLIGHT_FRAME - 25,
                HIGHLIGHT_FRAME + 25,
                color=highlight_color,
                alpha=0.3,
                zorder=0
            )

        # Utiliser les mêmes limites d'axe Y que le graphique global
        ax.set_ylim(global_y_limits)

        # Inversion de l'axe Y
        ax.invert_yaxis()

        # Ajouter des lignes horizontales pointillées
        ax.axhline(y=-1, color='red', linestyle='--', linewidth=1, alpha=0.7)
        ax.axhline(y=0.78, color='green', linestyle='--', linewidth=1, alpha=0.7)

        # Configuration du graphique
        ax.set_title(f'Évolution du contact {pair}\npour le résidu {residue_id}')
        ax.set_xlabel('Numéro de frame')
        ax.set_ylabel('Indice de frustration (inversé)')
        ax.grid(True)

        # Amélioration de l'affichage des ticks sur l'axe X
        ax.set_xticks(np.arange(0, 4901, 500))

        plt.tight_layout()

        # Sauvegarder le graphique
        pair_name = pair.replace(':', '_').replace('-', '_')
        plt.savefig(
            os.path.join(individual_dir, f"{os.path.basename(output_prefix)}_{pair_name}_individual.png"),
            bbox_inches='tight',
            dpi=300
        )
        plt.close()


def get_isolation_mode(dataframes_dir):
    """Détermine le mode d'isolation à partir du chemin du répertoire"""
    path_parts = os.path.normpath(dataframes_dir).split(os.sep)

    if 'True_isolated' in path_parts:
        return 'True_isolated'
    elif 'Isolated' in path_parts:
        return 'Isolated'
    else:
        return 'Not_isolated'


def main():
    parser = argparse.ArgumentParser(description='Analyse des données de frustration des contacts.')
    parser.add_argument('dataframes_dir', type=str,
                        help='Répertoire contenant les dataframes de contacts')
    parser.add_argument('--residue', type=str,
                        help='Résidu spécifique à analyser (format: Chain:AAResNum)')
    parser.add_argument('--frames', nargs='+', type=int,
                        help='Frames spécifiques à analyser (ex: 0 10)')
    parser.add_argument('--only-intra', action='store_true',
                        help='Afficher seulement les contacts intra-chaîne')
    parser.add_argument('--only-inter', action='store_true',
                        help='Afficher seulement les contacts inter-chaînes')

    args = parser.parse_args()

    # Vérification des arguments
    if args.only_intra and args.only_inter:
        print("Error: Vous ne pouvez pas utiliser --only-intra et --only-inter simultanément")
        return

    # Déterminer le type de contact à afficher
    contact_type = None
    if args.only_intra:
        contact_type = 'intra'
    elif args.only_inter:
        contact_type = 'inter'

    # Charge les données
    contact_data = load_contact_data(args.dataframes_dir)
    if contact_data is None:
        return

    # Extraire le nom de la protéine et le mode d'isolation
    protein_name = os.path.basename(os.path.normpath(args.dataframes_dir))
    isolation_mode = get_isolation_mode(args.dataframes_dir)

    # Crée le répertoire de sortie avec la structure correcte
    output_dir = os.path.join('../contact_analysis', isolation_mode, protein_name)
    os.makedirs(output_dir, exist_ok=True)

    # Analyse un résidu spécifique si demandé
    if args.residue:
        # Filtrer les contacts selon le type demandé
        filtered_contacts = filter_contacts(contact_data, args.residue, contact_type)

        # Analyser les contacts filtrés
        contacts, stats = analyze_residue_contacts(filtered_contacts, args.residue)

        if contacts is not None:
            print(f"\nAnalyse pour le résidu {args.residue}:")
            print(f"Contacts totaux: {stats['total_contacts']}")
            print(f"  - Intra-chaîne: {stats['intra_chain']}")
            print(f"  - Inter-chaînes: {stats['inter_chain']}")
            print(f"Frustration moyenne: {stats['avg_frustration']:.3f}")
            print(f"Contacts très frustrés: {stats['highly_frustrated']}")
            print(f"Contacts peu frustrés: {stats['minimally_frustrated']}")
            print(f"Contacts neutres: {stats['neutral']}")
            print(f"Frames analysées: {stats['frames_analyzed']} (de {stats['first_frame']} à {stats['last_frame']})")

            # Sauvegarde les contacts dans un fichier
            residue_file = os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}.csv")
            contacts.to_csv(residue_file, index=False)
            print(f"\nContacts détaillés sauvegardés dans {residue_file}")

            # Trace l'évolution
            y_limits, color_map = plot_residue_contacts(
                contacts, args.residue,
                os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}"), args.dataframes_dir , contact_type ,
                frames=args.frames if args.frames else None
            )
            print(f"Graphique d'évolution sauvegardé")

            # Si on a tracé toutes les frames, génère aussi les graphiques individuels
            if args.frames is None and y_limits is not None:
                plot_individual_contacts(
                    contacts, args.residue,
                    os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}"),
                    color_map, y_limits
                )
                print("Graphiques individuels des contacts sauvegardés dans le dossier 'individual_contacts'")
    else:
        print("Veuillez spécifier un résidu avec --residue pour l'analyse")


if __name__ == "__main__":
    main()
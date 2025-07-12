import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from collections import defaultdict
import numpy as np

# Configuration des frames par défaut à visualiser
DEFAULT_FRAMES = [0, 4900]  # Modifier cette liste pour changer les frames par défaut


def load_contact_data(protein_name, isolate, true_isolate):
    """Charge les données de contacts précédemment sauvegardées"""

     # Determine the subdirectory based on isolation mode
    if true_isolate:
        subdir = "True_isolated"
    elif isolate:
        subdir = "Isolated"
    else:
        subdir = "Not_isolated"

    data_dir = os.path.join('../contact_dataframes', subdir, protein_name)

    if not os.path.exists(data_dir):
        print(f"Error: No data found for protein {protein_name}")
        return None

    contact_files = [f for f in os.listdir(data_dir) if f.endswith('.tsv')]
    if not contact_files:
        print(f"Error: No TSV files found in {data_dir}")
        return None

    # Charger tous les fichiers et les combiner
    dfs = []
    for f in contact_files:
        df = pd.read_csv(os.path.join(data_dir, f), sep='\t')
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


def plot_residue_contacts(contacts, residue_id, output_prefix, frames=None):
    """Visualisation des contacts pour un résidu spécifique"""
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
    plt.title(f'Évolution des contacts pour le résidu {residue_id}\nFrames: {frame_label}')
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

    # Ajouter la légende des tailles
    #plt.gca().add_artist(first_legend)
    #plt.legend(handles=size_legend,
    #           title='Type de contact',
    #           bbox_to_anchor=(1.05, 0.7),
    #           loc='upper left')

    plt.tight_layout()

    # Obtenir les limites de l'axe Y avant de sauvegarder
    y_limits = plt.gca().get_ylim()

    plt.savefig(f"{output_prefix}_contact_evolution_{frame_label}.png",
                bbox_inches='tight',
                dpi=300)
    plt.close()

    return y_limits, color_map


def plot_individual_contacts(contacts, residue_id, output_prefix, color_map, global_y_limits):
    """Génère des graphiques individuels pour chaque contact avec mise en évidence d'un frame spécifique"""
    # Frame à mettre en évidence (configurable facilement)
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
        #ax.invert_yaxis()

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

def main():
    parser = argparse.ArgumentParser(description='Analyse des données de frustration des contacts.')
    parser.add_argument('protein', type=str, help='Nom de la protéine à analyser')
    parser.add_argument('--residue', type=str, help='Résidu spécifique à analyser (format: Chain:AAResNum)')
    parser.add_argument('--frames', nargs='+', type=int, help='Frames spécifiques à analyser (ex: 0 10)')
    parser.add_argument('--only-intra', action='store_true', help='Afficher seulement les contacts intra-chaîne')
    parser.add_argument('--only-inter', action='store_true', help='Afficher seulement les contacts inter-chaînes')
    parser.add_argument('--isolated', action='store_true', help='Whether the frustration was calculated by separating the chains')
    parser.add_argument('--true_isolate', action='store_true', default=False,
                        help='Similar to isolate but results are saved in True_isolate directory')

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

    # Charge les données selon le mode spécifié
    contact_data = load_contact_data(args.protein, args.isolated, args.true_isolate)
    print(contact_data)
    if contact_data is None:
        return

     # Determine the subdirectory based on isolation mode
    if args.true_isolate:
        subdir = "True_isolated"
    elif args.isolate:
        subdir = "Isolated"
    else:
        subdir = "Not_isolated"


    # Crée le répertoire de sortie approprié
    output_dir = os.path.join('../contact_analysis', subdir, args.protein)
    os.makedirs(output_dir, exist_ok=True)

    # Analyse un résidu spécifique si demandé
    if args.residue:
        # Filtrer les contacts selon le type demandé
        filtered_contacts = filter_contacts(contact_data, args.residue, contact_type)

        # Analyser les contacts filtrés
        contacts, stats = analyze_residue_contacts(filtered_contacts, args.residue)

        if contacts is not None:
            print(f"\nAnalyse pour le résidu {args.residue}:")
            print(f"Mode: {'Isolated' if args.isolated else 'Not isolated'}")
            print(f"Contacts totaux: {stats['total_contacts']}")
            print(f"  - Intra-chaîne: {stats['intra_chain']}")
            print(f"  - Inter-chaînes: {stats['inter_chain']}")
            print(f"Frustration moyenne: {stats['avg_frustration']:.3f}")
            print(f"Contacts très frustrés: {stats['highly_frustrated']}")
            print(f"Contacts peu frustrés: {stats['minimally_frustrated']}")
            print(f"Contacts neutres: {stats['neutral']}")
            print(f"Frames analysées: {stats['frames_analyzed']} (de {stats['first_frame']} à {stats['last_frame']})")

            # Sauvegarde les contacts dans un fichier
            residue_file = os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}.tsv")
            contacts.to_csv(residue_file, sep='\t', index=False)
            print(f"\nContacts détaillés sauvegardés dans {residue_file}")

            # Trace l'évolution
            y_limits, color_map = plot_residue_contacts(
                contacts, args.residue,
                os.path.join(output_dir, f"residue_{args.residue.replace(':', '_')}"),
                frames=args.frames if args.frames else None
            )
            print(f"Graphique d'évolution sauvegardé")
            # Si on a tracé toutes les frames, génère aussi les graphiques individuels
            #y_limits = (np.float64(0.9734499999999999), np.float64(-2.05045))
            #color_map = {'0:?187-0:Y189': (0.9677975592919913, 0.44127456009157356, 0.5358103155058701), '0:Y189-0:I191': (0.8836443049112893, 0.5240073524369634, 0.19569304285113343), '0:Y189-0:E193': (0.710130687316902, 0.6046852192663268, 0.19426060163712158), '0:Y189-0:?194': (0.5432776721247529, 0.6540981095185215, 0.19324494273892204), '0:Y189-0:R197': (0.19592059105779686, 0.6981620017487838, 0.3452219818913641), '0:T184-0:Y189': (0.2067117296964458, 0.6829103404254792, 0.5829988925822328), '0:E183-0:Y189': (0.21420912437215422, 0.6714963557258681, 0.6986206664203177), '0:P146-0:Y189': (0.22537170008202412, 0.6531400148480775, 0.841007805313343), '0:Y189-0:R192': (0.5596943802099308, 0.5764402169887779, 0.9583930713150347), '0:P150-0:Y189': (0.8578978803740231, 0.44058452715322166, 0.957819659566579), '0:R147-0:Y189': (0.9628653850704806, 0.4025928454059796, 0.7779310354076443)}
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
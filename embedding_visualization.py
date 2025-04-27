import marimo

__generated_with = "0.10.2"
app = marimo.App(
    width="medium",
    app_title="Rhea_database_visualization",
    auto_download=["html"],
)


@app.cell
def _():
    #imports

    import pandas as pd
    import numpy as np
    np.random.seed(77)
    import more_itertools
    import os
    import io
    import re
    import h5py
    import gzip
    import requests
    import PIL
    import itertools
    import random
    import urllib.request
    import urllib.error
    import matplotlib.pyplot as plt
    import matplotlib.font_manager
    from matplotlib.colors import ListedColormap
    from matplotlib.cm import viridis

    from collections import Counter
    from fontTools import ttLib

    from datetime import datetime
    current_date = datetime.now().date()
    date = current_date.strftime('%Y-%m-%d')

    import umap

    import ipywidgets as widgets
    from IPython.display import display

    # from graphein.protein.utils import download_alphafold_structure

    import py3Dmol

    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, pairwise_distances

    from nltk import ngrams, FreqDist

    from scipy.spatial.distance import cdist
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests

    from tqdm import tqdm

    import marimo as mo

    # import altair as alt

    import rdkit

    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw

    import altair as alt
    from rdkit.Chem import rdmolfiles
    return (
        AllChem,
        Chem,
        Counter,
        Draw,
        FreqDist,
        KMeans,
        ListedColormap,
        PIL,
        alt,
        cdist,
        current_date,
        date,
        datetime,
        display,
        gzip,
        h5py,
        hypergeom,
        io,
        itertools,
        matplotlib,
        mo,
        more_itertools,
        multipletests,
        ngrams,
        np,
        os,
        pairwise_distances,
        pd,
        plt,
        py3Dmol,
        random,
        rdkit,
        rdmolfiles,
        re,
        requests,
        silhouette_score,
        tqdm,
        ttLib,
        umap,
        urllib,
        viridis,
        widgets,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##1. Read Data CSV for UMAP values""")
    return


@app.cell
def _(pd):
    data_df = pd.read_csv("umap_smiles.csv")
    data_df = data_df[['Compounds', 'SMILES', 'UMAP_1', 'UMAP_2']]
    data_df
    return (data_df,)


@app.cell
def _(KMeans, data_df):
    #Clustering - Number of clusters
    kmeans = KMeans(n_clusters=5, random_state=42)
    clusters = kmeans.fit_predict(data_df[['UMAP_1', 'UMAP_2']])
    data_df['clusters'] = clusters
    return clusters, kmeans


@app.cell
def _(alt):
    def scatter(df, pan = False):
        selection = alt.selection_interval(bind='scales', translate=pan, zoom=True)
        return (alt.Chart(df)
            .mark_circle()
            .encode(
                x=alt.X("UMAP_1:Q"),
                y=alt.Y("UMAP_2:Q"),
                color=alt.Color("clusters:N")
            )
            .properties(width=1000, height=500)).add_params(selection)
    return (scatter,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##2. Data Visualization of Rhea database Steroids""")
    return


@app.cell(hide_code=True)
def _(mo):
    checkbox = mo.ui.checkbox(label="Toggle Pan")
    checkbox
    return (checkbox,)


@app.cell
def _(checkbox, data_df, mo, scatter):
    chart = mo.ui.altair_chart(scatter(data_df, checkbox.value))
    chart
    return (chart,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""##3. Data values from visualization""")
    return


@app.cell(hide_code=True)
def _(Chem, Draw, display, table):
    for smile_val in table.value["SMILES"]:
        # Check if there is a semicolon in the SMILES string
        if ";" in smile_val:
            # Split the SMILES into individual components
            smile_parts = smile_val.split(";")
            # Convert each part into a molecule
            molecules = [Chem.MolFromSmiles(smile) for smile in smile_parts]
            # Display molecules side by side
            img = Draw.MolsToGridImage(molecules, molsPerRow=len(molecules))
            display(img)
        else:
            # If no semicolon, convert the SMILES and display as normal
            molecules = Chem.MolFromSmiles(smile_val)
            img = Draw.MolToImage(molecules)
            display(img)
    return img, molecules, smile_parts, smile_val


@app.cell(hide_code=True)
def _(chart, mo):
    table = mo.ui.table(chart.value)
    table
    return (table,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""##4. 3d Visualization (Up to 9 structures)""")
    return


@app.cell
def _(mo):
    button = mo.ui.run_button(label = "Download HTML")
    button
    return (button,)


@app.cell
def _(AllChem, Chem, button, os, py3Dmol, rdmolfiles, table):
    # Example list of SMILES strings (replace with your chart.value['SMILES'])
    smiles_list = table.value['SMILES']

    # Keep only the first 9 (py3Dmol grid max is 3x3)
    selected_smiles = smiles_list[:9]

    if button.value:
        # Step 1: Generate PDBs from SMILES
        pdb_list = []
        for smi in selected_smiles:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            pdb = rdmolfiles.MolToPDBBlock(mol)
            pdb_list.append((smi, pdb))

        # Step 2: Set up py3Dmol viewer with 3x3 grid
        viewer = py3Dmol.view(viewergrid=(3,3), width=900, height=900)
        viewer.setBackgroundColor('white')

        for i, (smi, pdb) in enumerate(pdb_list):
            row = i // 3
            col = i % 3
            viewer.addModel(pdb, 'pdb', viewer=(row, col))
            viewer.setStyle({'stick': {}}, viewer=(row, col))
            viewer.zoomTo(viewer=(row, col))
            viewer.addLabel(smi, {'position': {'x': 0, 'y': -5, 'z': 0}}, viewer=(row, col))

        # Step 3: Generate HTML
        html_code = f"""
        <html>
          <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol.js"></script>
          </head>
          <body>
            {viewer._make_html()}
          </body>
        </html>
        """

        # Step 4: Save to unique filename if file already exists
        base_filename = "3D_rhea_values"
        extension = ".html"
        filename = base_filename + extension
        counter = 1
        while os.path.exists(filename):
            filename = f"{base_filename}_{counter}{extension}"
            counter += 1

        with open(filename, "w") as f:
            f.write(html_code)

        print(f"✅ HTML with 3D structures saved as '{filename}'")
    return (
        base_filename,
        col,
        counter,
        extension,
        f,
        filename,
        html_code,
        i,
        mol,
        pdb,
        pdb_list,
        row,
        selected_smiles,
        smi,
        smiles_list,
        viewer,
    )


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

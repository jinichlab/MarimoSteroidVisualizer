import marimo

__generated_with = "0.10.2"
app = marimo.App()


@app.cell
def simple_ui():
    import marimo as mo
    import pandas as pd
    from sklearn.cluster import KMeans
    import altair as alt
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdmolfiles
    from IPython.display import display, HTML
    import base64
    import py3Dmol
    return (
        AllChem,
        Chem,
        Draw,
        HTML,
        KMeans,
        alt,
        base64,
        display,
        mo,
        pd,
        py3Dmol,
        rdmolfiles,
    )


@app.cell
def _(pd):
    data_df = pd.read_csv("umap_smiles.csv")
    data_df = data_df[['Compounds', 'SMILES', 'UMAP_1', 'UMAP_2']]
    # data_df
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
            ).add_params(selection)
    return (scatter,)


@app.cell
def _(mo):
    checkbox = mo.ui.checkbox(label="Toggle Pan")
    checkbox
    return (checkbox,)


@app.cell
def _(checkbox, data_df, mo, scatter):
    chart = mo.ui.altair_chart(scatter(data_df, checkbox.value))
    chart
    return (chart,)


@app.cell
def _(Chem, Draw, HTML, display, pd, table):
    for _, rows in table.value.iterrows():
        smile_val = rows["SMILES"]
        compound_name = rows["Compounds"]
        width = 1000
        height = 1000

        if pd.isna(smile_val):
            continue

        if ";" in smile_val:
            smile_parts = smile_val.split(";")
            molecules = [Chem.MolFromSmiles(smile) for smile in smile_parts]
            img = Draw.MolsToGridImage(molecules, molsPerRow=len(molecules), subImgSize=(width, height))
        else:
            molecule = Chem.MolFromSmiles(smile_val)
            img = Draw.MolToImage(molecule, size=(width, height))

        # Display the compound name and the image
        display(HTML(f"<p><b>{compound_name}</b></p>"))
        display(img)

    return (
        compound_name,
        height,
        img,
        molecule,
        molecules,
        rows,
        smile_parts,
        smile_val,
        width,
    )


@app.cell
def _():
    # table.value.index
    return


@app.cell
def _(chart, mo):
    table = mo.ui.table(chart.value)
    table
    return (table,)


@app.cell
def _(mo):
    button = mo.ui.run_button(label = "Generate 3D Structures")
    button
    return (button,)


@app.cell
def _(AllChem, Chem, base64, button, py3Dmol, rdmolfiles, table):
    smiles_list = table.value['SMILES']

    # Keep only the first 9 (py3Dmol grid max is 3x3)
    selected_smiles = smiles_list[:9]
    download_link = ''
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

        b64_html = base64.b64encode(html_code.encode()).decode()
        
        # Step 3: Create the download link
        download_link = f'<a download="3D_rhea_values.html" href="data:text/html;base64,{b64_html}">Download 3D Viewer HTML</a>'
    return (
        b64_html,
        col,
        download_link,
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
def _(download_link, mo):
    mo.md(download_link)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

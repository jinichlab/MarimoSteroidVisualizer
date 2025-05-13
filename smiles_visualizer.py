import marimo

__generated_with = "0.10.2"
app = marimo.App(width="medium")


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
    raw_data_df = pd.read_csv("umap_compound_smiles.csv")
    # raw_data_df
    return (raw_data_df,)


@app.cell
def _(pd, raw_data_df):
    # chebi_df = pd.read_csv("compounds.tsv", sep="\t", encoding="ISO-8859-1", usecols=["ID", "NAME"])
    chebi_df = pd.read_csv("chebi_lookup_minimal.csv")
    chebi_df["ID"] = chebi_df["ID"].astype(str)
    raw_data_df["ChEBI ID"] = raw_data_df["ChEBI ID"].astype(str)

    data_df = raw_data_df.merge(chebi_df, left_on="ChEBI ID", right_on="ID", how="left")
    data_df = data_df.rename(columns={"NAME": "Compound"}).drop(columns=["ID"])
    data_df = data_df.drop("SMILES_clean", axis=1)
    return chebi_df, data_df


@app.cell
def _(KMeans, data_df):
    #Clustering - Number of clusters
    kmeans = KMeans(n_clusters=15, random_state=42)
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
def _():
    # table.value.index
    return


@app.cell
def _(chart, mo):
    table = mo.ui.table(chart.value)
    table
    return (table,)


@app.cell
def _(Chem, Draw, HTML, chebi_df, display, pd, table):
    # Setup: Lookup dictionary for ChEBI names
    chebi_df["ID"] = chebi_df["ID"].astype(str).str.strip()
    chebi_lookup = dict(zip(chebi_df["ID"], chebi_df["NAME"]))

    # Font and image settings
    width, height = 400, 400  # smaller individual molecule image size
    mols_per_row = 2          # better readability

    # Loop through each row in your table
    for _, row in table.value.iterrows():
        smile_val = row["SMILES"]
        chebi_val = row["ChEBI ID"]
        protein_name = row.get("Protein names", "[No Protein Name]")
        entry_name = row.get("Entry Name", "[No Entry Name]")
        organism_name = row.get("Organism", "[No Organism Name]")

        if pd.isna(smile_val) or pd.isna(chebi_val):
            continue

        # Display protein name above molecule(s)
        display(HTML(f"<h3 style='font-size:23px;'><strong>Protein Name:</strong> {protein_name}</h3>"))
        display(HTML(f"<h3 style='font-size:23px;'><strong>Entry Name:</strong> {entry_name}</h3>"))
        display(HTML(f"<h3 style='font-size:23px;'><strong>Organism:</strong> {organism_name}</h3>"))




        if ";" in smile_val:
            smile_parts = [s.strip() for s in smile_val.split(";")]
            chebi_parts = [c.strip() for c in chebi_val.split(";")]

            # Convert SMILES → RDKit Mol
            molecules = [Chem.MolFromSmiles(smi) for smi in smile_parts]

            # Get names using ChEBI lookup
            compound_names = [chebi_lookup.get(cid, f"[Unknown: {cid}]") for cid in chebi_parts]

            # Draw molecules with legends
            img = Draw.MolsToGridImage(
                molecules,
                molsPerRow=mols_per_row,
                subImgSize=(width, height),
                legends=compound_names,
                legendFontSize=24
            )
            display(img)
            display(HTML("<hr style='border:1px solid #ccc;'>"))

        else:
            molecule = Chem.MolFromSmiles(smile_val)
            chebi_id = str(chebi_val).strip()
            compound_name = chebi_lookup.get(chebi_id, f"[Unknown: {chebi_id}]")

            img = Draw.MolToImage(molecule, size=(width, height), legend=compound_name)
            display(img)
            display(HTML("<hr style='border:1px solid #ccc;'>"))
    return (
        chebi_id,
        chebi_lookup,
        chebi_parts,
        chebi_val,
        compound_name,
        compound_names,
        entry_name,
        height,
        img,
        molecule,
        molecules,
        mols_per_row,
        organism_name,
        protein_name,
        row,
        smile_parts,
        smile_val,
        width,
    )


@app.cell
def _():
    return


@app.cell
def _(mo):
    button = mo.ui.run_button(label = "Generate 3D Structures")
    button
    return (button,)


@app.cell
def _(
    AllChem,
    Chem,
    base64,
    button,
    chebi_df,
    pd,
    py3Dmol,
    rdmolfiles,
    table,
):
    # 1. Build ChEBI ID → Name dictionary (safe var name)
    chebi_name_lookup = dict(zip(chebi_df["ID"].astype(str).str.strip(), chebi_df["NAME"]))

    # 2. Get first 9 SMILES and ChEBI ID entries
    mol3d_smiles_list = table.value['SMILES']
    mol3d_chebi_ids = table.value['ChEBI ID']
    mol3d_selected_smiles = mol3d_smiles_list[:9]
    mol3d_selected_chebis = mol3d_chebi_ids[:9]

    mol3d_download_link = ''

    if button.value:
        mol3d_pdb_data = []

        for mol3d_smi_str, mol3d_chebi_str in zip(mol3d_selected_smiles, mol3d_selected_chebis):
            if pd.isna(mol3d_smi_str) or pd.isna(mol3d_chebi_str):
                continue

            mol3d_smi_parts = [s.strip() for s in mol3d_smi_str.split(";")]
            mol3d_chebi_parts = [c.strip() for c in mol3d_chebi_str.split(";")]

            for smi, cid in zip(mol3d_smi_parts, mol3d_chebi_parts):
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    print(f"⚠️ Invalid SMILES skipped: {smi}")
                    continue
                mol = Chem.AddHs(mol)
                try:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)
                    pdb = rdmolfiles.MolToPDBBlock(mol)
                    label = chebi_name_lookup.get(cid, f"[Unknown: {cid}]")
                    mol3d_pdb_data.append((label, pdb))
                except Exception as e:
                    print(f"⚠️ Error processing {smi}:\n{e}")

        # 3. Create py3Dmol viewer (max 3x3 grid)
        mol3d_total = min(len(mol3d_pdb_data), 9)
        mol3d_rows = (mol3d_total + 2) // 3
        mol3d_viewer = py3Dmol.view(viewergrid=(mol3d_rows, 3), width=900, height=900)
        mol3d_viewer.setBackgroundColor('white')

        for i, (label, pdb) in enumerate(mol3d_pdb_data[:9]):
            row_idx = i // 3
            col_idx = i % 3
            mol3d_viewer.addModel(pdb, 'pdb', viewer=(row_idx, col_idx))
            mol3d_viewer.setStyle({'stick': {}}, viewer=(row_idx, col_idx))
            mol3d_viewer.zoomTo(viewer=(row_idx, col_idx))
            mol3d_viewer.addLabel(label, {
                'position': {'x': 0, 'y': 5, 'z': 0},
                'fontSize': 16,
                'fontColor': 'black',
                'backgroundOpacity': 0,  # transparent background
                'inFront': True
            }, viewer=(row_idx, col_idx))

        # 4. Generate downloadable HTML
        mol3d_html_code = f"""
        <html>
          <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol.js"></script>
          </head>
          <body>
            {mol3d_viewer._make_html()}
          </body>
        </html>
        """
        mol3d_b64_html = base64.b64encode(mol3d_html_code.encode()).decode()
        mol3d_download_link = f'<a download="3D_rhea_compounds.html" href="data:text/html;base64,{mol3d_b64_html}">Download 3D Viewer HTML</a>'
    return (
        chebi_name_lookup,
        cid,
        col_idx,
        i,
        label,
        mol,
        mol3d_b64_html,
        mol3d_chebi_ids,
        mol3d_chebi_parts,
        mol3d_chebi_str,
        mol3d_download_link,
        mol3d_html_code,
        mol3d_pdb_data,
        mol3d_rows,
        mol3d_selected_chebis,
        mol3d_selected_smiles,
        mol3d_smi_parts,
        mol3d_smi_str,
        mol3d_smiles_list,
        mol3d_total,
        mol3d_viewer,
        pdb,
        row_idx,
        smi,
    )


@app.cell
def _(mo, mol3d_download_link):
    mo.md(mol3d_download_link)
    return


@app.cell
def _(chebi_df, data_df):
    # 1. Extract all ChEBI IDs from your dataset
    all_chebi_ids = data_df["ChEBI ID"].dropna().astype(str)

    # 2. Flatten multi-ID rows like "16113; 46898"
    chebi_id_list = []
    for entry in all_chebi_ids:
        chebi_id_list.extend([cid.strip() for cid in entry.split(";")])

    # 3. Deduplicate
    unique_chebi_ids = sorted(set(chebi_id_list))

    # 4. Look up names from the full chebi_df
    chebi_df["ID"] = chebi_df["ID"].astype(str).str.strip()
    chebi_subset = chebi_df[chebi_df["ID"].isin(unique_chebi_ids)][["ID", "NAME"]].copy()

    # 5. Write to CSV
    chebi_subset.to_csv("chebi_lookup_minimal.csv", index=False)
    return (
        all_chebi_ids,
        chebi_id_list,
        chebi_subset,
        entry,
        unique_chebi_ids,
    )


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

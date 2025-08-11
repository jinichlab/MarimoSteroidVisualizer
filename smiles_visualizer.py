import marimo

__generated_with = "0.14.12"
app = marimo.App(width="medium")


@app.cell
def simple_ui():
    import marimo as mo
    import pandas as pd
    from sklearn.cluster import KMeans
    import altair as alt
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, rdmolfiles
    from IPython.display import display, HTML
    import base64
    import py3Dmol

    # For LLM summarization and secure API usage
    from dotenv import load_dotenv
    import os
    import openai
    import instructor
    from anthropic import Anthropic
    from pydantic import BaseModel
    from typing import List
    from jinja2 import Template
    import re
    import pyarrow

    from io import BytesIO

    import ast
    return (
        AllChem,
        BaseModel,
        BytesIO,
        Chem,
        Draw,
        HTML,
        KMeans,
        List,
        alt,
        ast,
        base64,
        display,
        load_dotenv,
        mo,
        openai,
        os,
        pd,
        py3Dmol,
        rdmolfiles,
        re,
    )


@app.cell
def _():
    import requests
    def get_uniprot_entry(uniprot_id):
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception("Failed to fetch UniProt data")
        return response.json()

    # Example usage
    ent = get_uniprot_entry("P12345")
    print(ent["organism"]["scientificName"])
    print(ent["sequence"]["value"])
    return


@app.cell
def _(pd):
    small_molecule_df = pd.read_csv("small_molecule.csv")
    # raw_data_df
    return (small_molecule_df,)


@app.cell
def _(pd):
    # protein_embedding_df = pd.read_csv("protein_embeddings.csv")
    protein_embedding_df = pd.read_csv("sequence_embeddings.csv")
    protein_embedding_df.drop('embedding', axis=1, inplace=True)
    return (protein_embedding_df,)


@app.cell
def _(display):
    display('Select the type of graph')
    return


@app.cell
def _(mo):
    dropdown = mo.ui.dropdown(["small molecule centric", "protein centric"])
    dropdown
    return (dropdown,)


@app.cell
def _(pd):
    full_chebi_df = pd.read_csv("compounds.tsv", sep="\t", encoding="ISO-8859-1", usecols=["ID", "NAME"])
    return (full_chebi_df,)


@app.cell
def _(full_chebi_df):
    full_chebi_df['ID'] = full_chebi_df['ID'].apply(lambda x: int(x))
    return


@app.cell
def _(ast, full_chebi_df, small_molecule_df):
    # Select relevant columns
    small_molecule_names = small_molecule_df[['SMILES', 'ChEBI ID']].copy()

    # Safely extract the number from the string list
    small_molecule_names['ChEBI ID'] = small_molecule_names['ChEBI ID'].apply(
        lambda x: int(ast.literal_eval(x)[0])
    )

    # Now merge with full_chebi_df
    merged_df = small_molecule_names.merge(
        full_chebi_df,
        left_on='ChEBI ID',
        right_on='ID',
        how='left'
    )
    return (merged_df,)


@app.cell
def _(merged_df):
    merged_df[['SMILES', 'ChEBI ID', 'NAME']].to_csv('small_molecules_names.csv')
    return


@app.cell
def _(ast, dropdown, pd, protein_embedding_df, small_molecule_df):
    # chebi_df = pd.read_csv("compounds.tsv", sep="\t", encoding="ISO-8859-1", usecols=["ID", "NAME"])
    chebi_df = pd.read_csv("chebi_lookup_minimal.csv")
    chebi_df["ID"] = chebi_df["ID"].astype(int)

    raw_data_df = pd.DataFrame()
    if dropdown.value == "small molecule centric":
        raw_data_df = small_molecule_df.copy()
        raw_data_df['ChEBI ID'] = raw_data_df['ChEBI ID'].apply(lambda x: int(ast.literal_eval(x)[0]))
        data_df = raw_data_df.merge(chebi_df, left_on="ChEBI ID", right_on="ID", how="left")
        data_df = data_df.rename(columns={"NAME": "Compound"}).drop(columns=["ID"])
    elif dropdown.value == "protein centric":
        data_df = protein_embedding_df.copy()
    return chebi_df, data_df


@app.class_definition
class SkipCell(Exception):
    pass


@app.cell
def _():
    # if dropdown.value != None:
    #     display(data_df)
    return


@app.cell
def _(KMeans, data_df, dropdown):
    #Clustering - Number of clusters
    if dropdown.value !=None:
        kmeans = KMeans(n_clusters=5, random_state=42)
        clusters = kmeans.fit_predict(data_df[['UMAP_1', 'UMAP_2']])
        data_df['clusters'] = clusters
    return


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
            ).add_params(selection).properties(width = 800, height = 600)
    return (scatter,)


@app.cell
def _(mo):
    checkbox = mo.ui.checkbox(label="Toggle Pan")
    return (checkbox,)


@app.cell
def _(checkbox, data_df, display, dropdown, mo, scatter):
    if dropdown.value != None:
        chart = mo.ui.altair_chart(scatter(data_df, checkbox.value))
        display(chart)
    return (chart,)


@app.cell
def _(checkbox):
    checkbox
    return


@app.cell
def _(chart, display, dropdown, mo):
    if dropdown.value != None:
        table = mo.ui.table(chart.value)
        display(table)
    return (table,)


@app.cell(hide_code=True)
def rdkittohtml(BytesIO, Draw, base64, chebi_df):
    # Helper to convert RDKit image to base64 HTML

    # Setup: Lookup dictionary for ChEBI names
    chebi_df["ID"] = chebi_df["ID"].astype(str).str.strip()
    chebi_lookup = dict(zip(chebi_df["ID"], chebi_df["NAME"]))

    # Font and image settings
    width, height = 400, 400  # individual molecule image size

    max_scroll_height = 200   # scroll box height
    mols_per_row = 2          # readability
    def mols_to_base64_html(molecules, legends=None, mols_per_row=mols_per_row):
        img = Draw.MolsToGridImage(
            molecules,
            molsPerRow=mols_per_row,
            subImgSize=(width, height),
            legends=legends,
            legendFontSize=24
        )
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        b64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
        return f'<img src="data:image/png;base64,{b64}" />'
    return chebi_lookup, max_scroll_height, mols_to_base64_html


@app.function(hide_code=True)
def counter(protein_list):
    items = []
    for i, p in enumerate(set(protein_list), start=1):
        # give each <li> some bottom margin for a blank line effect
        items.append(
          f'<li style="margin-bottom:15px;">'
          f'<b>{i}</b>. {p}'
          '</li>'
        )
    return "".join(items)


@app.cell(hide_code=True)
def _(
    Chem,
    HTML,
    ast,
    chebi_lookup,
    display,
    dropdown,
    max_scroll_height,
    mols_to_base64_html,
):
    def display_compound_with_scroll(smile_val, chebi_val, proteins, entries):
        """
        Render one SMILES entry as:
         • An RDKit-drawn grid of molecule images (with legends from ChEBI)
         • A fixed-height, scrollable box listing the associated proteins with Entries
        """
        # 1) Normalize & dedupe proteins + entries
        if len(proteins) == 1 and isinstance(proteins[0], str) and ";" in proteins[0]:
            proteins = [p.strip() for p in proteins[0].split(";")]
        if len(entries) == 1 and isinstance(entries[0], str) and ";" in entries[0]:
            entries = [e.strip() for e in entries[0].split(";")]

        # Deduplicate by (protein, entry) pair
        seen = set()
        paired = []
        for prot, entry in zip(proteins, entries):
            key = (prot.strip(), entry.strip())
            if key not in seen:
                seen.add(key)
                paired.append(key)

        # 2) Prepare molecule images
        smile_parts = [s.strip() for s in str(smile_val).split(";")]
        mols = [Chem.MolFromSmiles(s) for s in smile_parts]
        chebi_ids = [c.strip() for c in str(chebi_val).split(";")]
        legends = [chebi_lookup.get(cid, f"[Unknown:{cid}]") for cid in chebi_ids]
        img_html = mols_to_base64_html(mols, legends)

        # 3) Build the numbered <ul>
        if dropdown.value == 'small molecule centric' and isinstance(paired[0][0], str):
            paired = [ast.literal_eval(p) if isinstance(p, str) else p for p in paired]

        protein_items = ""
        for i, (prot, entry) in enumerate(paired, 1):
            protein_items += f"""
            <li style="margin-bottom: 12px;">
                <strong>{i}.</strong> <strong>Entry:</strong> {entry}<br> <strong>Name:</strong> {prot} <br><a href="https://alphafold.ebi.ac.uk/entry/{entry}"
               target="_blank"
               style="color: #1a73e8; text-decoration: underline;">
               AlphaFold Structure
            </a>
            """
        protein_list_html = f"<ul>{protein_items}</ul>"

        # 4) Wrap in a scroll box
        scroll_box = f"""
          <div style="
             max-height: {max_scroll_height}px;
             overflow-y: auto;
             border: 1px solid #ccc;
             padding: 8px;
             width: 350px;
             background: #fafafa;
             word-wrap: break-word;
          ">
            {protein_list_html}
          </div>
        """

        # 5) Compose the final layout
        html = f"""
        <div style="
          display: flex;
          gap: 20px;
          align-items: flex-start;
          margin-bottom: 24px;
        ">
          <div>{img_html}</div>
          <div>
            <div style="font-size:16px; margin-bottom:4px;">
              <strong>Proteins:</strong> (scroll)
            </div>
            {scroll_box}
          </div>
        </div>
        """
        display(HTML(html))

    return (display_compound_with_scroll,)


@app.cell
def _(HTML, ast, display, display_compound_with_scroll, dropdown, pd, table):

    if dropdown.value != None:
        # Loop through each row in your table
        for _, row in table.value.iterrows():
            smile_val = row["SMILES"]
            chebi_val = row["ChEBI ID"]

            proteins = row.get("Protein names", [])
            entries = row.get("Entry", [])

            # Handle semicolon-delimited or stringified lists
            if dropdown.value == "small molecule centric":
                if isinstance(proteins, str):
                    proteins = ast.literal_eval(proteins)
                if isinstance(entries, str):
                    entries = ast.literal_eval(entries)
            else:
                if isinstance(proteins, str):
                    proteins = [p.strip() for p in proteins.split(";")]
                if isinstance(entries, str):
                    entries = [e.strip() for e in entries.split(";")]

            if pd.isna(smile_val) or pd.isna(chebi_val):
                continue

            display_compound_with_scroll(smile_val, chebi_val, proteins, entries)
            display(HTML('<hr style="border:1px solid #ccc; margin:16px 0;">'))

    return


@app.cell
def _(display, dropdown, mo):
    if dropdown.value != None:
        if len(dropdown.value)>0:
            button = mo.ui.run_button(label = "Generate Small molecule 3D Structures")
        display(button)
    return (button,)


@app.cell
def _(
    AllChem,
    Chem,
    base64,
    button,
    chebi_df,
    dropdown,
    pd,
    py3Dmol,
    rdmolfiles,
    table,
):
    if dropdown.value != None:
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
    return (mol3d_download_link,)


@app.cell
def _(dropdown, mo, mol3d_download_link):
    if dropdown.value != None:
        mo.md(mol3d_download_link)
    return


@app.cell
def _(chebi_df, data_df, dropdown):
    if dropdown.value != None:
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
    return


@app.cell
def _():
    return


@app.cell
def _(BaseModel, List):
    class OutputFormat(BaseModel):
        summary: str
        highlights: List[str]
        tldr: str
    return


@app.cell
def _(mo):
    text_input = mo.ui.text(label="Enter Prompt Here:")

    mo.md(f"""
    AI generated information:

    {text_input} 

    """).batch(text_input=text_input).form()
    return (text_input,)


@app.cell
def _(load_dotenv, os):
    env_path = os.path.join(os.path.dirname(__file__), ".env")
    load_dotenv(dotenv_path=env_path)
    'TOGETHER_API_KEY' in str(os.environ)
    TOGETHER_KEY = os.getenv("TOGETHER_API_KEY")
    return (TOGETHER_KEY,)


@app.cell
def _(TOGETHER_KEY, mo, openai, text_input):
    client = openai.OpenAI(
        api_key=TOGETHER_KEY,
        base_url="https://api.together.xyz/v1"
    )

    if text_input.value.strip():
        with mo.status.spinner("Generating information..."):
            response = client.chat.completions.create(
                model="mistralai/Mixtral-8x7B-Instruct-v0.1",
                messages=[
                    {
                        "role": "user",
                        "content": f"""Answer the prompt below.

    Your answer must begin with 'TL;DR:' on a new line followed by a one-sentence summary. Then provide 3–5 bullet point highlights.

    Prompt:
    {text_input.value}
    """
                    }
                ],
                max_tokens=1000
            )

            output_text = response.choices[0].message.content
            rendered = f"""## Response from Mixtral\n{output_text}"""
    else:
        output_text = ""
        rendered = "⚠️ No input provided."

    # output_text  # ← Add this so the next cell can access it
    return (output_text,)


@app.cell
def _(mo, output_text, re):
    match = re.search(
        r"^TL;DR[:：]?\s*(.+?)(?:\n{2,}|\n(?=[\d\-•]))",  # captures TLDR before list
        output_text.strip(),
        re.DOTALL | re.IGNORECASE
    )

    if match:
        summary = match.group(1).strip()
        main_content = output_text[match.end():].strip()
    else:
        summary = ""
        main_content = output_text.strip()

    # Only render if there is content
    formatted_output = f"""
    ## TL;DR
    {summary}

    ## Key Highlights
    {main_content}
    """ if len(main_content) + len(summary) > 0 else ""

    mo.md(formatted_output)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

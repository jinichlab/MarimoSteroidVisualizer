import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium", app_title="Steroid Visualizer")


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
    import openai
    from pydantic import BaseModel
    from typing import List
    from io import BytesIO
    import ast
    from openai import OpenAI
    import os, json, numpy as np, faiss

    return (
        AllChem,
        BaseModel,
        BytesIO,
        Chem,
        Draw,
        HTML,
        KMeans,
        List,
        OpenAI,
        alt,
        ast,
        base64,
        display,
        faiss,
        json,
        load_dotenv,
        mo,
        np,
        openai,
        os,
        pd,
        py3Dmol,
        rdmolfiles,
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
    protein_embedding_df = protein_embedding_df.drop(columns = ['Unnamed: 0'])
    # protein_embedding_df.drop('embedding', axis=1, inplace=True)
    return (protein_embedding_df,)


@app.cell
def _(mo):
    enter_button = mo.ui.run_button(label = "Enter App", full_width=True)
    return (enter_button,)


@app.cell(hide_code=True)
def _(mo):
    intro_text = """
    # 🧭 Welcome to the Steroid–Protein Navigator

    This portal lets you explore Nature’s steroids and the proteins and enzymes they interact with—all through an interactive interface.

    ## Views available

    **• Steroid-centric view →** Explore the chemical space of steroids.  
    Each point represents a steroid molecule; clicking one highlights the interacting proteins.

    **• Protein-centric view →** Explore the protein landscape.  
    Each point represents a steroid-binding protein; selecting one reveals the steroids or sterol-like molecules it associates with.

    ## How to use it

    You can navigate either view by selecting clusters or individual entries.  
    When you click on a point, a detailed table allowing you to select values and view structures, names, and alphafold structures.

    A built-in language model companion (trained on 100+ DeepResearch reports)
    helps explain the biology, biochemistry, and known steroid–protein interactions
    associated with each cluster.

    Use the Navigator to uncover patterns, explore relationships, and generate new ideas.

    <span style="font-size: 0.9em; font-style: italic; color: #555;">
    Tip: To search for specific values, choose a view and select all clusters. Then use the table’s search bar to look up specific entries.
    </span>
    """
    intro_screen = mo.md(intro_text)
    intro_screen = mo.vstack([intro_screen], align="center")

    return (intro_screen,)


@app.cell(hide_code=True)
def _(mo):
    loading_spinner = mo.Html("""
    <div id="spinner-wrapper" style="
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        height: 200px;
        gap: 12px;
    ">

      <!-- Circle Loader -->
      <div style="
          width: 48px;
          height: 48px;
          border: 6px solid #e0e0e0;
          border-top-color: #0077ff;
          border-radius: 50%;
          animation: spin 0.8s linear infinite;
      "></div>

      <div style="font-size: 16px; color: #444;">
        Loading…
      </div>

    </div>

    <style>
    @keyframes spin {
      to { transform: rotate(360deg); }
    }
    </style>

    <script>
    setTimeout(() => {
        const el = document.getElementById("spinner-wrapper");
        if (el) el.style.display = "none";
    }, 2000);   // <-- X seconds (2000 ms)
    </script>
    """)

    return


@app.cell
def _(display, enter_button, intro_screen):
    if not enter_button.value:
        display(intro_screen)
        display(enter_button)
    else:
        None

    return


@app.cell
def _(mo):
    temp = mo.ui.run_button(label="Enter App")


    temp = temp

    return


@app.cell
def _(display, enter_button, mo):
    dropdown = mo.ui.dropdown(["protein centric", "small molecule centric"])
    if enter_button.value:
        display('Select the type of graph')
        display(dropdown)  
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
        data_df = data_df.rename(columns={"NAME": "Compound Name"}).drop(columns=["ID"])
    elif dropdown.value == "protein centric":
        data_df = protein_embedding_df
    return chebi_df, data_df


@app.cell
def _(KMeans, data_df, dropdown):
    #Clustering - Number of clusters
    if dropdown.value !=None:
        kmeans = KMeans(n_clusters=9, random_state=42)
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
def _(checkbox, enter_button):
    if enter_button.value:
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
    return max_scroll_height, mols_to_base64_html


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
def _(Chem, HTML, display, max_scroll_height, mols_to_base64_html):
    def display_compound_with_scroll(smile_val, compound_name_val, proteins, entries):
        """Show stacked (Name → SMILES image) pairs on the left and a scrollable protein list on the right."""

        # Proteins / entries: split if they came as single semicolon-joined strings, then dedupe pairs
        if len(proteins) == 1 and isinstance(proteins[0], str) and ";" in proteins[0]:
            proteins = [p.strip() for p in proteins[0].split(";")]
        if len(entries) == 1 and isinstance(entries[0], str) and ";" in entries[0]:
            entries = [e.strip() for e in entries[0].split(";")]

        seen, paired = set(), []
        for prot, entry in zip(proteins, entries):
            key = (prot.strip(), entry.strip())
            if key not in seen:
                seen.add(key)
                paired.append(key)

        # Parse SMILES / names and align lengths
        smile_parts = [s.strip() for s in str(smile_val).split(";") if s and str(s).strip()]
        name_parts  = [n.strip() for n in str(compound_name_val).split(";") if n and str(n).strip()]
        mols = [Chem.MolFromSmiles(s) for s in smile_parts]

        if len(name_parts) < len(mols):
            name_parts += ["[Name missing]"] * (len(mols) - len(name_parts))
        else:
            name_parts = name_parts[:len(mols)]

        # Build Name : Image blocks 
        pair_blocks = []
        for nm, m in zip(name_parts, mols):
            img_html = (
                '<div style="font-size:12px;color:#999;">[Invalid SMILES]</div>'
                if m is None else mols_to_base64_html([m], [""])
            )
            pair_blocks.append(
                f"""
                <div style="margin-bottom:14px;">
                  <div style="font-size:15px; margin-bottom:6px;"><strong>Compound Name:</strong> {nm.title()}</div>
                  <div>{img_html}</div>
                </div>
                """
            )
        small_molecule_html = "".join(pair_blocks)

        # Protein list (numbered, scrollable)
        protein_items = "".join(
            f"""
            <li style="margin-bottom:10px;">
              <strong>{i}.</strong> <strong>Entry:</strong> {e}<br>
              <strong>Name:</strong> {p}<br>
              <a href="https://alphafold.ebi.ac.uk/entry/{e}" target="_blank" style="text-decoration:underline; color:#1a73e8;">
                AlphaFold Structure
              </a>
            </li>
            """
            for i, (p, e) in enumerate(paired, 1)
        )
        protein_box_html = f"""
          <div style="font-size:16px; margin-bottom:6px;"><strong>Proteins</strong></div>
          <div style="max-height:{max_scroll_height}px; overflow-y:auto; border:1px solid #e0e0e0; padding:8px; background:#fafafa;">
            <ul style="margin:0; padding-left:18px;">{protein_items}</ul>
          </div>
        """

        # Layout
        html = f"""
        <div style="display:grid; grid-template-columns:1fr 360px; gap:20px; align-items:start; margin-bottom:24px;">
          <div>{small_molecule_html}</div>
          <div>{protein_box_html}</div>
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
            compound_name = row["Compound Name"]

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

            if pd.isna(smile_val) or pd.isna(compound_name):
                continue

            display_compound_with_scroll(smile_val, compound_name, proteins, entries)
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
    #3D Small molecule download

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


@app.cell(hide_code=True)
def _(BaseModel, List):
    class OutputFormat(BaseModel):
        summary: str
        highlights: List[str]
        tldr: str
    return


@app.cell
def _(load_dotenv, openai, os):
    #Create openai instance, load API key

    load_dotenv(dotenv_path=".env")
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
    if not OPENAI_API_KEY:
        raise ValueError("❌ OPENAI_API_KEY not found in .env file")

    # Create OpenAI client

    client = openai.Client(api_key=OPENAI_API_KEY)
    return (client,)


@app.cell
def _(OpenAI, faiss, json, np, os):
    #Load RAG data and embed texts

    # retriever.py — load index & catalog once; simple top-k search

    STORE = "rag_store"
    _index = faiss.read_index(os.path.join(STORE, "index.faiss"))
    _catalog = [json.loads(l) for l in open(os.path.join(STORE, "catalog.jsonl"), "r", encoding="utf-8")]
    _client = OpenAI()
    INDEX_PATH = os.path.join(STORE, "index.faiss")
    CATALOG_PATH = os.path.join(STORE, "catalog.jsonl")
    # Load FAISS index once
    index = faiss.read_index(INDEX_PATH)

    def _embed(texts):
        resp = _client.embeddings.create(model="text-embedding-3-large", input=texts)
        arr = np.array([d.embedding for d in resp.data], dtype="float32")
        arr /= (np.linalg.norm(arr, axis=1, keepdims=True) + 1e-12)
        return arr
    return CATALOG_PATH, index


@app.cell
def _(CATALOG_PATH, client, index, json, np):
    #Find closest chunks and embed query in vector space

    # Load catalog (maps vector IDs → text/metadata)
    with open(CATALOG_PATH, "r", encoding="utf-8") as f:
        catalog = [json.loads(line) for line in f]

    def embed_query(query: str) -> np.ndarray:
        """Turn query into normalized vector"""
        resp = client.embeddings.create(
            model="text-embedding-3-large",
            input=[query]
        )
        v = np.array(resp.data[0].embedding, dtype="float32")
        return v / (np.linalg.norm(v) + 1e-12)

    def retrieve(query: str, k: int = 6):
        """Search index for top-k chunks matching query"""
        qv = embed_query(query).reshape(1, -1)
        scores, idxs = index.search(qv, k)
        hits = []
        for score, idx in zip(scores[0], idxs[0]):
            if idx == -1:
                continue
            row = catalog[idx]
            row["score"] = float(score)
            hits.append(row)
        return hits
    return (retrieve,)


@app.cell
def _(client, dropdown, np, retrieve, table):
    #Actual model with system message. Also finds similarity with RAG data

    def my_model(messages):
        # Keep history like before
        history = []
        last_role = None
        for m in messages:
            role = "assistant" if m.role == "assistant" else "user"
            if role == last_role:
                continue
            history.append({"role": role, "content": m.content})
            last_role = role
        if not history:
            history = [{"role": "user", "content": "Hello!"}]




        selection_keywords = (
            "select", "selected", "selection",
            "selected values", "values selected",
            "the selected", "from the table", "from the list",
            "highlighted", "chosen", "picked"
        )
        # Get the latest user query
        user_q = next((m["content"] for m in reversed(history) if m["role"]=="user"), "Hello!")
        selected = any(kw in user_q.lower() for kw in selection_keywords)

        system_msg = (
            "You are a precise assistant. "
            "When context is available, always prioritize it. "
            "If no context is relevant, answer briefly (2–3 simple sentences) from your own knowledge and note that it is general information. "
            "Always mention what document is being referenced. If no document is referenced mention using general information.")

        if dropdown.value != None:
            if selected:
                if len(np.array(table.value)) == 0:
                    system_msg = "If asked about selected proteins, respond by requesting the user to first select values from the table."
                elif dropdown.value == "small molecule centric":
                    user_q = f"{user_q}: {np.array(table.value['Compound'])}. "

                elif dropdown.value == "protein centric":
                    user_q = f"{user_q}: {np.array(table.value['Protein names'])}. "
        else:
            system_msg = "If asked about selected proteins, respond by requesting the user to first select type of graph from the dropdown."



        # RAG part (Retrieval):
        hits = retrieve(user_q, k=6)
        context_blocks = []
        for i, h in enumerate(hits, 1):
            header = f"[{i}] {h['paper']} — {h.get('section') or '(no section)'} — page {h.get('page')}"
            preview = h["text"]
            context_blocks.append(f"{header}\n{preview}")
        context = "\n\n".join(context_blocks) if context_blocks else "(no relevant context found)"

        # Call the chat completion
        resp = client.chat.completions.create(
            model="gpt-4o-mini",  # or gpt-4.1-mini if you want
            messages=[
                {"role": "system", "content": system_msg},
                {"role": "user", "content": f"Context:\n{context}\n\nQuestion: {user_q}"}
            ],
            temperature=0.2,
        )

        return resp.choices[0].message.content

    return (my_model,)


@app.cell
def _(display, enter_button, mo, my_model):
    chat_ui = mo.ui.chat(
        my_model,
        prompts=[
            "Summarize this chat app in one sentence.",
            "Tell me about the values selected",
            "Tell me about steroid-enzyme interaction in general"
        ],
    )
    # Only the messages area scrolls
    chat_ui.style({
        "height": "400px",
        "overflow-y": "auto",
    })
    if enter_button.value:
        display(mo.Html(
            f"""
            <div id="chat-container" style="height:400px; overflow-y:auto; border:1px solid #ccc; border-radius:8px; padding:8px;">
                {chat_ui}
            </div>
            <script>
            const container = document.getElementById("chat-container");
            const observer = new MutationObserver(() => {{
                container.scrollTop = container.scrollHeight;
            }});
            observer.observe(container, {{ childList: true, subtree: true }});
            </script>
            """
        ))


    return


@app.cell
def _(display, enter_button, mo):
    if enter_button.value:
        display(mo.md(
            """
    
        <span style='font-size:16px'><em>Hint: Unsure what to ask? Click the prompt icon next to the input bar for ideas<em></span>
        """
        ))
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

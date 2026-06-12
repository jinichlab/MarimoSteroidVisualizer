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
    import difflib
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
        difflib,
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
    small_molecule_df = pd.read_csv("small_molecule_centric.csv")
    cols = small_molecule_df.columns.tolist()
    cols = ["Compound Name"] + [c for c in cols if c != "Compound Name"]
    small_molecule_df = small_molecule_df[cols]
    # small_molecule_df.drop(columns="Unnamed: 0")
    # raw_data_df
    return (small_molecule_df,)


@app.cell
def _(pd):
    # protein_embedding_df = pd.read_csv("protein_embeddings.csv")
    protein_embedding_df = pd.read_csv(
        "protein_sequence_embedding.csv",
        dtype={"ChEBI ID": str, "Rhea ID": str, "SMILES": str},
    )
    protein_embedding_df = protein_embedding_df.drop(columns = ['Unnamed: 0'])
    # NaN -> "" for ID/SMILES columns so .split() and dict lookups behave
    for _c in ("ChEBI ID", "Rhea ID", "SMILES"):
        if _c in protein_embedding_df.columns:
            protein_embedding_df[_c] = protein_embedding_df[_c].fillna("")
    # protein_embedding_df.drop('embedding', axis=1, inplace=True)
    return (protein_embedding_df,)


@app.cell
def _(mo):
    enter_button = mo.ui.run_button(label = "Enter App", full_width=True)
    help_button = mo.ui.run_button(label = "Help")
    return enter_button, help_button


@app.cell
def _(mo):
    b1 = mo.ui.run_button(label = 'b1')
    b2 = mo.ui.run_button(label = 'b2')
    return


@app.cell(hide_code=True)
def _(mo):
    intro_text = mo.Html("""
    <div style="text-align:center; margin-top:10px;">

    <span style="font-size:34px; font-weight:600;">
    🧭 Welcome to Nature's Steroid Atlas
    </span>

    <p style="max-width:750px; margin:12px auto 0; font-size:15px;">
    Explore Nature’s steroids and the proteins they interact with — all through an interactive, unified interface.
    </p>

    <hr style="width:60%; margin:16px auto;">

    <h3 style="margin:0;">Views Available</h3>

    <div style="text-align:left; max-width:650px; margin:10px auto; font-size:15px;">
    <b>• Steroid-centric view →</b> Explore the chemical space of steroids.<br>
    <span style="margin-left:18px;">Each point represents a steroid molecule; clicking one highlights interacting proteins.</span>
    <br><br>
    <b>• Protein-centric view →</b> Explore the protein landscape.<br>
    <span style="margin-left:18px;">Each point represents a steroid-binding protein; selecting one reveals associated steroids.</span>
    <br><br>
    <b>• Natural and Synthetic steroid view →</b> Explore known natural steroids and synthetic steroids
    <br>
    <span style="margin-left:18px;">
    Each point represents a steroid; selecting one reveals associated proteins when available. 
    <br>
    <span style="margin-left:18px;">
    Use this view to explore synthetic steroids by comparing them to known natural steroids.  
    </span>


    <h3 style="margin-top:12px;">How to Use It</h3>

    <p style="text-align:left; max-width:650px; margin:10px auto; font-size:15px;">
    You can navigate either view by selecting clusters or clicking on individual points to explore the underlying steroids or proteins. 
    When you select a point, a detailed table appears where you can then select specific rows to examine structures, compound or protein names, associated identifiers 
    (such as UniProt, ChEBI, or Rhea IDs), and AlphaFold protein models. A built-in language model companion—trained on more than 130 
    DeepResearch reports-helps explain the biology, biochemistry, and known steroid–protein interactions connected to the selected entries. 
    Use the Navigator to uncover patterns, investigate molecular relationships, and generate new ideas as you explore.
    </p>

    </div>
    """)

    intro_screen = mo.vstack([intro_text], align="center")

    tip = mo.Html("""<p style="font-size:13px; opacity:0.7;">
    Tip: To search for specific values, choose a view and select all clusters. Then use the table’s search bar to look up specific entries.
    </p>""")

    return (intro_screen,)


@app.cell
def _(display, enter_button, help_button, intro_screen):
    if (not enter_button.value) or help_button.value:
        display(intro_screen)
        display(enter_button)
    else:
        None
    return


@app.cell
def _(display, enter_button, help_button, mo):
    dropdown = mo.ui.dropdown(["protein centric", "small molecule centric", "Natural and Synthetic steroids"])
    if enter_button.value:
        display(help_button)
        display('Select the type of graph')
        display(dropdown)  
    return (dropdown,)


@app.cell
def _(display, dropdown, mo, steroids):
    #Search Specific <protein/small molecule> — by name OR UniProt entry
    # Default so downstream cells can always read query.value without
    # NameError on the landing page (before any dropdown selection).
    query = mo.ui.text(placeholder="")
    if dropdown.value!=None:
        if dropdown.value == "protein centric":
            placeholder = "enter protein name, gene, or UniProt entry (e.g. P12345)"
        elif dropdown.value == "small molecule centric":
            placeholder = "enter compound name or UniProt entry of a binding protein"
        else:
            placeholder = f"enter {steroids} or UniProt entry of a binding protein"
        query = mo.ui.text(placeholder=placeholder)
        display(query)
    return (query,)


@app.cell
def _(data_df, difflib, display, dropdown, query):
    if dropdown.value!=None and len(query.value)>0:
        # Search across name + UniProt-accession columns. For molecule views
        # the Entry column is a list-repr string like "['P12345','Q67890']";
        # substring match still picks up the accession inside.
        if dropdown.value == 'protein centric':
            primary_col = "Protein names"
            search_cols = ["Protein names", "Entry", "Entry Name", "Gene Names"]
        else:
            primary_col = "Compound Name"
            search_cols = ["Compound Name", "Entry"]
        search_cols = [c for c in search_cols if c in data_df.columns]

        query_str = query.value.strip().lower()

        # 1. substring match across every search column (OR), case-insensitive
        mask = data_df[search_cols[0]].astype(str).str.lower().str.contains(
            query_str, na=False)
        for _c in search_cols[1:]:
            mask = mask | data_df[_c].astype(str).str.lower().str.contains(
                query_str, na=False)
        exact_matches = data_df[mask]

        if len(exact_matches) > 0:
            filtered_df = exact_matches
        else:
            # 2. fuzzy fallback on the primary NAME column only — accessions
            #    shouldn't be fuzzy-matched (P12345 vs P12346 is nonsense).
            choices = data_df[primary_col].astype(str).unique()
            nearest = difflib.get_close_matches(query.value, choices,
                                                n=5, cutoff=0.3)
            if len(nearest) > 0:
                filtered_df = data_df[data_df[primary_col].astype(str).isin(nearest)]
            else:
                filtered_df = data_df.iloc[0:0]

        display("Closest matches to entry:")
        # Show name + Entry side-by-side so the user can confirm what matched.
        # Surface is_new so a hit on a newly-recruited protein is obvious.
        _show_cols = [c for c in [primary_col, "Entry", "Gene Names",
                                   "clusters", "is_new"]
                      if c in filtered_df.columns]
        display(filtered_df[_show_cols])
    return


@app.cell
def _():
    #Names

    # full_chebi_df = pd.read_csv("Names ref/chebi_compound_names.tsv", sep="\t", encoding="ISO-8859-1", usecols=["ID", "NAME"])
    # full_chebi_df['ID'] = full_chebi_df['ID'].apply(lambda x: int(x))
    return


@app.cell
def _():
    # Cleaning small molecule + Names

    # # Select relevant columns
    # small_molecule_names = small_molecule_df[['SMILES', 'ChEBI ID']].copy()

    # # Safely extract the number from the string list
    # small_molecule_names['ChEBI ID'] = small_molecule_names['ChEBI ID'].apply(
    #     lambda x: int(ast.literal_eval(x)[0])
    # )

    # # Now merge with full_chebi_df
    # merged_df = small_molecule_names.merge(
    #     full_chebi_df,
    #     left_on='ChEBI ID',
    #     right_on='ID',
    #     how='left'
    # )

    # merged_df[['SMILES', 'ChEBI ID', 'NAME']].to_csv('small_molecules_names.csv')
    return


@app.cell
def _():
    # Combined df cleanup

    # combined_df = pd.read_csv("natural_synthetic_steroids.csv")
    # combined_df = combined_df.rename(columns = {'type': 'clusters', 'Name':'Compound Name'})
    # combined_df = combined_df.drop(columns=['Unnamed: 0'])
    # list_cols = ["Protein names", "Entry", "Entry Name", "Gene Names", "ChEBI ID"]

    # def fix_list_cell(x):
    #     if pd.isna(x):
    #         return []                   # NaN → empty list
    #     if isinstance(x, list):
    #         return x                    # already good
    #     if isinstance(x, str):
    #         try:
    #             return ast.literal_eval(x)  # "['A','B']" → ['A','B']
    #         except:
    #             return [x]              # fallback: wrap single string
    #     return [x]                      # fallback for anything else

    # for col in list_cols:
    #     if col in combined_df.columns:
    #         combined_df[col] = combined_df[col].apply(fix_list_cell)

    return


@app.cell
def _(dropdown, pd, protein_embedding_df, small_molecule_df):
    # chebi_df = pd.read_csv("compounds.tsv", sep="\t", encoding="ISO-8859-1", usecols=["ID", "NAME"])
    chebi_df = pd.read_csv("Names ref/chebi_lookup_minimal.csv")
    chebi_df["ID"] = chebi_df["ID"].astype(int)

    raw_data_df = pd.DataFrame()
    # Default so downstream cells never see NameError on the landing page
    # (before any dropdown selection). All downstream cells already gate
    # rendering on `dropdown.value != None`, so an empty frame is safe.
    data_df = pd.DataFrame()
    if dropdown.value == "small molecule centric":
        data_df = small_molecule_df
        data_df = data_df.drop(columns="Unnamed: 0")
    elif dropdown.value == "protein centric":
        data_df = protein_embedding_df
    elif dropdown.value == "Natural and Synthetic steroids":
        data_df = pd.read_csv('natural_synthetic_steroids.csv')
    return chebi_df, data_df


@app.cell
def _(KMeans, data_df, dropdown):
    # Clusters are precomputed in the build scripts: KMeans on the 2D UMAP
    # coords with k chosen automatically by silhouette score (not a fixed 9).
    # Natural/Synthetic keeps its class label in `clusters`, so skip it.
    if dropdown.value not in (None, "Natural and Synthetic steroids"):
        if "clusters" not in data_df.columns or data_df["clusters"].isna().all():
            # defensive fallback only if the precomputed column is missing
            data_df["clusters"] = KMeans(
                n_clusters=9, random_state=42
            ).fit_predict(data_df[["UMAP_1", "UMAP_2"]])
    return


@app.cell
def _(alt):
    import colorsys

    def _build_palette(n):
        """Generate n maximally-spread distinct hex colours.

        Strategy: golden-ratio hue stepping (avoids the clumping you get
        from a regular hue sweep) combined with 3 alternating lightness
        levels, so adjacent indices vary in BOTH hue and lightness.
        Gives ~25-30 perceptually distinct colours; beyond that some
        pairs will still look similar — that is a perceptual limit of
        the human eye, not a palette problem.
        """
        phi_inv = 0.6180339887498949   # golden-ratio conjugate
        lightnesses = [0.55, 0.40, 0.72]
        sat = 0.70
        out = []
        for i in range(max(n, 1)):
            h = (i * phi_inv) % 1.0
            l = lightnesses[i % len(lightnesses)]
            r, g, b = colorsys.hls_to_rgb(h, l, sat)
            out.append(f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}")
        return out

    def scatter(df, pan=False):
        selection = alt.selection_interval(bind='scales', translate=pan, zoom=True)

        # Use Altair's qualitative scheme for small k (looks cleaner with
        # familiar Tableau colours), and an explicit hand-built palette
        # for many clusters so every category gets a UNIQUE hex code.
        try:
            n = int(df["clusters"].nunique())
        except Exception:
            n = 10
        if n <= 10:
            scale = alt.Scale(scheme="tableau10")
        elif n <= 20:
            scale = alt.Scale(scheme="tableau20")
        else:
            domain = sorted(df["clusters"].unique().tolist(), key=lambda x: str(x))
            scale = alt.Scale(domain=domain, range=_build_palette(len(domain)))

        # Always show the legend. For many clusters, lay it out in multiple
        # short columns with a small symbol so it fits next to the plot.
        if n <= 20:
            legend = alt.Legend(title="cluster")
        else:
            legend = alt.Legend(
                title="cluster",
                columns=2,
                symbolLimit=200,
                labelFontSize=9,
                titleFontSize=10,
                symbolSize=60,
            )
        tooltip = [c for c in ("Entry", "Protein names", "Compound Name",
                               "Gene Names", "Annotation", "Paper", "clusters")
                   if c in df.columns]

        # 5-pointed star SVG path for newly-recruited proteins so they
        # visually pop against the dense circle background.
        STAR = ("M0,-1 L0.22,-0.31 L0.95,-0.31 L0.36,0.12 L0.59,0.81 "
                "L0,0.4 L-0.59,0.81 L-0.36,0.12 L-0.95,-0.31 L-0.22,-0.31 Z")

        has_new = "is_new" in df.columns and (df["is_new"] == 1).any()

        # If `is_match` column is present and some rows are 0, render the
        # search highlight using *conditional encoding inside a single
        # chart* (rather than layering two charts). marimo's altair_chart
        # wrapper injects its own selection signal at a fixed name and
        # that signal does not resolve on a layered spec.
        search_active = ("is_match" in df.columns
                         and (df["is_match"] == 0).any())

        new_pred = "datum.is_new == 1"

        if not search_active:
            # Base case: every point colored by cluster.
            # New proteins drawn as larger black-edged stars, ordered last.
            if has_new:
                shape_enc = alt.condition(new_pred, alt.value(STAR), alt.value("circle"))
                size_enc  = alt.condition(new_pred, alt.value(220), alt.value(18))
                stroke_enc = alt.condition(new_pred, alt.value("black"), alt.value("white"))
                stroke_w_enc = alt.condition(new_pred, alt.value(1.2), alt.value(0.15))
                order_enc = alt.Order("is_new:Q", sort="ascending")
            else:
                shape_enc = alt.value("circle")
                size_enc = alt.value(18)
                stroke_enc = alt.value("white")
                stroke_w_enc = alt.value(0.15)
                order_enc = alt.Undefined

            return (
                alt.Chart(df)
                .mark_point(filled=True, opacity=0.85)
                .encode(
                    x=alt.X("UMAP_1:Q"), y=alt.Y("UMAP_2:Q"),
                    color=alt.Color("clusters:N", scale=scale, legend=legend),
                    shape=shape_enc,
                    size=size_enc,
                    stroke=stroke_enc,
                    strokeWidth=stroke_w_enc,
                    order=order_enc,
                    tooltip=tooltip,
                )
                .add_params(selection)
                .properties(width=800, height=600)
            )

        # Search active: dim non-matches; matches get color + larger size.
        # New proteins still keep their star shape, whether matched or not.
        match_pred = "datum.is_match == 1"
        if has_new:
            shape_enc = alt.condition(new_pred, alt.value(STAR), alt.value("circle"))
            # Stars are larger than circles in both states
            size_enc = alt.condition(
                f"({match_pred}) && ({new_pred})", alt.value(320),
                alt.condition(match_pred, alt.value(160),
                              alt.condition(new_pred, alt.value(120), alt.value(8)))
            )
        else:
            shape_enc = alt.value("circle")
            size_enc = alt.condition(match_pred, alt.value(160), alt.value(8))

        return (
            alt.Chart(df)
            .mark_point(filled=True, stroke="black", strokeWidth=0.6)
            .encode(
                x=alt.X("UMAP_1:Q"), y=alt.Y("UMAP_2:Q"),
                color=alt.condition(
                    match_pred,
                    alt.Color("clusters:N", scale=scale, legend=legend),
                    alt.value("#cccccc"),
                ),
                shape=shape_enc,
                size=size_enc,
                opacity=alt.condition(match_pred, alt.value(1.0), alt.value(0.22)),
                # Draw matches last (on top) so they aren't covered by dim points
                order=alt.Order("is_match:Q", sort="ascending"),
                tooltip=tooltip,
            )
            .add_params(selection)
            .properties(width=800, height=600)
        )
    return (scatter,)


@app.cell
def _(mo):
    checkbox = mo.ui.checkbox(label="Toggle Pan")
    return (checkbox,)


@app.cell
def _(checkbox, data_df, display, dropdown, mo, query, scatter):
    # Build the df that gets plotted. If the search box has text, mark
    # matching rows with is_match=1 so scatter() can highlight them on
    # the UMAP (dim grey base + highlighted matches on top). When the
    # box is empty, render the chart normally.
    _df_to_plot = data_df.copy()
    _q = (query.value or "").strip().lower()
    if _q:
        if dropdown.value == "protein centric":
            _cols = ["Protein names", "Entry", "Entry Name", "Gene Names"]
        else:
            _cols = ["Compound Name", "Entry"]
        _cols = [_c for _c in _cols if _c in _df_to_plot.columns]
        if _cols:
            _mask = _df_to_plot[_cols[0]].astype(str).str.lower().str.contains(
                _q, na=False)
            for _c in _cols[1:]:
                _mask = _mask | _df_to_plot[_c].astype(str).str.lower().str.contains(
                    _q, na=False)
            _df_to_plot["is_match"] = _mask.astype(int)

    if dropdown.value != None:
        chart = mo.ui.altair_chart(scatter(_df_to_plot, checkbox.value))
        if _q and "is_match" in _df_to_plot.columns:
            _n_hits = int(_df_to_plot["is_match"].sum())
            display(f"Highlighting {_n_hits:,} match"
                    f"{'es' if _n_hits != 1 else ''} for "
                    f"'{query.value}' on the UMAP. "
                    "Drag to select clusters to explore them.")
        else:
            display("Drag and select clusters to explore them")
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
        # Render the Paper column as a clickable link if present.
        # format_mapping is keyed by column name and accepts a callable
        # that returns a marimo-renderable (mo.md works for markdown).
        _fmt = {}
        if "Paper" in chart.value.columns:
            def _fmt_paper(v):
                if not v or str(v).strip() == "":
                    return ""
                return mo.md(f"[Open paper]({v})")
            _fmt["Paper"] = _fmt_paper
        table = mo.ui.table(chart.value, format_mapping=_fmt) if _fmt else mo.ui.table(chart.value)
        display("Select values from the table to explore interacting proteins, and the structures")
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
def _(Chem, HTML, display, dropdown, max_scroll_height, mols_to_base64_html):
    def display_compound_with_scroll(smile_val, compound_name_val, proteins, entries, clusters):
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
          <div style="font-size:16px; margin-bottom:6px;">
            <strong>Proteins associated with this steroid ({len(paired)})</strong>
          </div>
          <div style="max-height:{max_scroll_height}px; overflow-y:auto; border:1px solid #e0e0e0; padding:8px; background:#fafafa;">
            <ul style="margin:0; padding-left:18px;">{protein_items}</ul>
          </div>
        """
        if dropdown.value == "Natural and Synthetic steroids" and clusters:
            cluster_html = f"""
            <div style="font-size:18px; font-weight:bold; margin-bottom:10px;">
              Type: {clusters.title()}
            </div>
            """
        else:
            cluster_html = ""

        # Layout
        html = f"""
        {cluster_html}
        <div style="display:grid; grid-template-columns:1fr 360px; gap:20px; align-items:start; margin-bottom:24px;">
          <div>{small_molecule_html}</div>
          <div>{protein_box_html}</div>
        </div>
        """
        display(HTML(html))
    return (display_compound_with_scroll,)


@app.cell
def _(Chem, HTML, display, mols_to_base64_html):
    def protein_view(row):
        """
        Protein-centric display for one row.
        Shows protein metadata + interacting small molecules (2 per row).
        """

        # ---------- Extract basic fields ----------
        proteins = row.get("Protein names", [])
        entries = row.get("Entry", [])
        smiles = row.get("SMILES", "")
        compound_names = row.get("Compound Name", "")
        organism = row.get("Organism", "[Unknown organism]")
        paper_url = str(row.get("Paper", "") or "").strip()
        annotation = str(row.get("Annotation", "") or "").strip()
        is_new = int(row.get("is_new", 0) or 0)

        # handle semicolon-separated or list values
        if isinstance(proteins, str):
            proteins = [p.strip() for p in proteins.split(";")]
        if isinstance(entries, str):
            entries = [e.strip() for e in entries.split(";")]

        # exactly one per row in protein-centric
        protein_name = proteins[0] if proteins else "[Unknown protein]"
        entry = entries[0] if entries else "[Unknown entry]"
        alphafold_link = f"https://alphafold.ebi.ac.uk/entry/{entry}"

        # ---------- Parse molecules ----------
        smile_parts = [s.strip() for s in str(smiles).split(";") if s.strip()]
        name_parts = [n.strip() for n in str(compound_names).split(";") if n.strip()]
        mols = [Chem.MolFromSmiles(s) for s in smile_parts]

        if len(name_parts) < len(mols):
            name_parts += ["[Name missing]"] * (len(mols) - len(name_parts))
        else:
            name_parts = name_parts[:len(mols)]

        # ---------- Build 2-per-row Molecule Cards ----------
        molecule_cards = []
        for nm, m in zip(name_parts, mols):

            if m is None:
                img_html = '<div style="font-size:12px;color:#999;">[Invalid SMILES]</div>'
            else:
                # smaller image
                img_html = mols_to_base64_html([m], [""], mols_per_row=1).replace("400", "250")

            molecule_cards.append(f"""
            <div style="width:48%; margin-bottom:16px; border:1px solid #ddd; padding:6px; border-radius:6px;">
                <div style="font-size:14px; font-weight:600; margin-bottom:4px;">
                    {nm.title()}
                </div>
                {img_html}
            </div>
            """)

        molecule_grid = f"""
        <div style="
            display:flex;
            flex-wrap:wrap;
            gap:10px;
            justify-content:space-between;
            margin-top:8px;">
            {''.join(molecule_cards)}
        </div>
        """

        # ---------- Final Layout ----------
        html = f"""
        <div style="padding:12px;">

            <div style="font-size:16px; margin-bottom:6px;">
                <strong>ENTRY:</strong> {entry}
            </div>

            <div style="font-size:16px; margin-bottom:6px;">
                <strong>Name:</strong> {protein_name}
            </div>

            <div style="font-size:16px; margin-bottom:6px;">
                <strong>Organism:</strong> {organism}
            </div>

            <div style="font-size:16px; margin-bottom:6px;">
                <strong>AlphaFold Structure:</strong>
                <a href="{alphafold_link}" target="_blank" style="color:#1a73e8; text-decoration:underline;">
                    {alphafold_link}
                </a>
            </div>

            {(
                '<div style="font-size:16px; margin-bottom:6px; padding:8px 10px; '
                'background:#fff8d6; border-left:4px solid #d4a017; border-radius:4px;">'
                '<strong>★ Newly recruited (2024-2026 literature)</strong><br>'
                f'<span style="font-size:13px; color:#444;">{annotation}</span>'
                '</div>'
            ) if is_new else ''}

            {(
                '<div style="font-size:16px; margin-bottom:20px;">'
                '<strong>Source paper:</strong> '
                f'<a href="{paper_url}" target="_blank" '
                'style="color:#1a73e8; text-decoration:underline;">'
                f'{paper_url}</a></div>'
            ) if paper_url else '<div style="margin-bottom:14px;"></div>'}

            <div style="font-size:18px; font-weight:bold; margin-bottom:10px;">
                Interacting Small Molecules:
            </div>

            {molecule_grid}

        </div>
        """

        display(HTML(html))

    return (protein_view,)


@app.cell
def _(
    HTML,
    ast,
    display,
    display_compound_with_scroll,
    dropdown,
    pd,
    protein_view,
    table,
):


    if dropdown.value != None:
        # Loop through each row in your table
        for _, row in table.value.iterrows():
            smile_val = row["SMILES"]
            compound_name = row["Compound Name"]

            proteins = row.get("Protein names", [])
            entries = row.get("Entry", [])
            cluster = row.get("clusters", [])  # 👈 get cluster here

            # Handle semicolon-delimited or stringified lists
            if dropdown.value == "protein centric":
                if isinstance(proteins, str):
                    proteins = [p.strip() for p in proteins.split(";")]
                if isinstance(entries, str):
                    entries = [e.strip() for e in entries.split(";")]
            else:
                if isinstance(proteins, str):
                    proteins = ast.literal_eval(proteins)
                if isinstance(entries, str):
                    entries = ast.literal_eval(entries)


            if pd.isna(smile_val) or pd.isna(compound_name):
                continue

            if dropdown.value == "protein centric":
                protein_view(row)   # 👈 new function
            else:
                display_compound_with_scroll(
                    smile_val, compound_name, proteins, entries, cluster
                )

            display(HTML('<hr style="border:1px solid #ccc; margin:16px 0;">'))

    return


@app.cell
def _(display, dropdown, mo, table):
    if dropdown.value != None and len(table.value)>0 and dropdown.value!="Natural and Synthetic steroids":
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
    if dropdown.value != None and len(table.value) and dropdown.value!="Natural and Synthetic steroids":
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
            mol3d_download_link = f'<a download="3D_rhea_compounds.html" href="data:text/html;base64,{mol3d_b64_html}">Link to Download 3D Viewer HTML</a>'
    return (mol3d_download_link,)


@app.cell
def _(display, dropdown, mo, mol3d_download_link, table):
    #3D Small molecule download

    if dropdown.value != None and len(table.value) and dropdown.value!="Natural and Synthetic steroids":
        display(mo.md(mol3d_download_link))
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
def _(mo):
    # Runtime RAG backend selector.
    #   local   - both embeddings and chat via local Ollama (no key needed)
    #   openai  - both via OpenAI cloud (needs OPENAI_API_KEY in .env)
    #   hybrid  - local nomic embeddings + OpenAI chat completions
    # The dropdown is reactive: flipping it re-loads the right FAISS index
    # without restarting the visualizer.
    backend_choice = mo.ui.dropdown(
        options={
            "Local — llama3.1:8b (offline, no key)": "local",
            "OpenAI — gpt-4o-mini (cloud, needs key)": "openai",
            "Hybrid — local retrieval + OpenAI chat": "hybrid",
        },
        value="Local — llama3.1:8b (offline, no key)",
        label="RAG backend",
    )
    return (backend_choice,)


@app.cell
def _(backend_choice, load_dotenv, openai, os):
    # Build the two clients the RAG pipeline needs: one for embeddings and
    # one for chat completions. In local/openai modes they're the same
    # object; in hybrid they diverge (Ollama for embeddings, OpenAI for chat).
    load_dotenv(dotenv_path=".env")
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
    _mode = backend_choice.value or "local"

    def _ollama_client():
        import urllib.request
        try:
            with urllib.request.urlopen(
                "http://127.0.0.1:11434/api/tags", timeout=2
            ) as _r:
                _r.read()
            return openai.Client(
                base_url="http://127.0.0.1:11434/v1",
                api_key="ollama",
            )
        except Exception as _e:
            print(
                f"⚠️  Local Ollama not reachable at 127.0.0.1:11434 ({_e}) — "
                "start it with `/home/adsiordia/ollama/bin/ollama serve &`."
            )
            return None

    def _openai_client():
        if not OPENAI_API_KEY:
            print("⚠️  OPENAI_API_KEY missing — add it to .env for OpenAI/Hybrid modes")
            return None
        return openai.Client(api_key=OPENAI_API_KEY)

    if _mode == "openai":
        embed_client = _openai_client()
        chat_client = embed_client
        print("✅ RAG backend: OpenAI cloud (chat + embeddings)")
    elif _mode == "hybrid":
        embed_client = _ollama_client()
        chat_client = _openai_client()
        print("✅ RAG backend: hybrid (local nomic embeddings + OpenAI chat)")
    else:
        embed_client = _ollama_client()
        chat_client = embed_client
        print("✅ RAG backend: local Ollama (chat + embeddings)")
    return chat_client, embed_client


@app.cell
def _(backend_choice, faiss, os):
    # FAISS index dimension is locked by the embedder, so the store must
    # match: nomic embeddings → rag_store_local (768-d); OpenAI embeddings
    # → rag_store (3072-d). Hybrid uses nomic, so it stays on the local store.
    _mode = backend_choice.value or "local"
    if _mode == "openai":
        STORE = "RAG Training/rag_store"
    else:
        STORE = "RAG Training/rag_store_local"

    INDEX_PATH = os.path.join(STORE, "index.faiss")
    CATALOG_PATH = os.path.join(STORE, "catalog.jsonl")
    index = faiss.read_index(INDEX_PATH)
    return CATALOG_PATH, index


@app.cell
def _(CATALOG_PATH, backend_choice, embed_client, index, json, np):
    # Load catalog (maps vector IDs → text/metadata)
    with open(CATALOG_PATH, "r", encoding="utf-8") as f:
        catalog = [json.loads(line) for line in f]

    _mode = backend_choice.value or "local"
    EMBED_MODEL = "text-embedding-3-large" if _mode == "openai" else "nomic-embed-text"

    def embed_query(query: str) -> np.ndarray:
        resp = embed_client.embeddings.create(model=EMBED_MODEL, input=[query])
        v = np.array(resp.data[0].embedding, dtype="float32")
        return v / (np.linalg.norm(v) + 1e-12)

    def retrieve(query: str, k: int = 6):
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
def _(backend_choice, chat_client, dropdown, np, retrieve, table):
    # Chat function used by mo.ui.chat. Builds a RAG-augmented prompt and
    # — when no retrieved chunk is on-topic — asks the model to fall back
    # on its own training knowledge with an explicit "[General knowledge]"
    # tag so the user can tell cited vs. uncited answers apart.
    def my_model(messages):
        if chat_client is None:
            _mode = backend_choice.value or "local"
            if _mode in ("openai", "hybrid"):
                return (
                    "🛈 OpenAI chat is unavailable — no OPENAI_API_KEY is "
                    "configured. Either add a key to `.env` and reload, or "
                    "switch the RAG backend dropdown above to **Local**."
                )
            return (
                "🛈 Local Ollama is not running on this server, so the chat "
                "assistant is offline. Every other feature of the app still "
                "works. Start Ollama with: "
                "`/home/adsiordia/ollama/bin/ollama serve &`."
            )

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
        user_q = next((m["content"] for m in reversed(history) if m["role"] == "user"), "Hello!")
        selected = any(kw in user_q.lower() for kw in selection_keywords)

        system_msg = (
            "You are a precise steroid-biology research assistant. You will "
            "be given a Context block of excerpts retrieved from a library "
            "of DeepResearch reports. Rules:\n"
            "1. If the Context contains material that answers the question, "
            "build your answer from it and cite each fact with the [n] tags "
            "that prefix the excerpts. Name the source paper.\n"
            "2. If the Context is empty, weakly relevant, or unrelated to "
            "the question, answer in 3-6 sentences from your own training "
            "knowledge and BEGIN your answer with the literal tag "
            "'[General knowledge — no relevant report found]' so the user "
            "knows it is uncited.\n"
            "3. Never fabricate a citation. If you didn't actually use a "
            "retrieved excerpt, don't cite it.\n"
        )

        if dropdown.value is not None:
            if selected:
                if len(np.array(table.value)) == 0:
                    system_msg = (
                        "If asked about selected proteins, respond by requesting "
                        "the user to first select values from the table."
                    )
                elif dropdown.value == "small molecule centric":
                    user_q = f"{user_q}: {np.array(table.value['Compound Name'])}. "
                elif dropdown.value == "protein centric":
                    user_q = f"{user_q}: {np.array(table.value['Protein names'])}. "
        else:
            system_msg = (
                "If asked about selected proteins, respond by requesting the user "
                "to first select type of graph from the dropdown."
            )

        hits = retrieve(user_q, k=6)
        context_blocks = []
        for i, h in enumerate(hits, 1):
            header = f"[{i}] {h['paper']} — {h.get('section') or '(no section)'} — page {h.get('page')}"
            context_blocks.append(f"{header}\n{h['text']}")
        context = "\n\n".join(context_blocks) if context_blocks else "(no relevant context found)"

        _mode = backend_choice.value or "local"
        CHAT_MODEL = "gpt-4o-mini" if _mode in ("openai", "hybrid") else "llama3.1:8b"

        resp = chat_client.chat.completions.create(
            model=CHAT_MODEL,
            messages=[
                {"role": "system", "content": system_msg},
                {"role": "user", "content": f"Context:\n{context}\n\nQuestion: {user_q}"},
            ],
            temperature=0.2,
        )
        return resp.choices[0].message.content

    return (my_model,)


@app.cell
def _(backend_choice, display, dropdown, mo, my_model):
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
    if dropdown.value !=None:
        display(mo.md("**RAG chat assistant**"))
        display(backend_choice)
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
def _(display, dropdown, mo):
    if dropdown.value !=None:
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

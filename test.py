import marimo

__generated_with = "0.10.2"
app = marimo.App()


@app.cell
def simple_ui():
    import marimo as mo
    import pandas as pd
    return mo, pd


@app.cell
def _(pd):
    data_df = pd.read_csv("umap_smiles.csv")
    data_df = data_df[['Compounds', 'SMILES', 'UMAP_1', 'UMAP_2']]
    data_df
    return (data_df,)


@app.cell
def _():
        
    return


if __name__ == "__main__":
    app.run()

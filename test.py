import marimo

__generated_with = "0.10.2"
app = marimo.App()


@app.cell
def simple_ui():
    import marimo as mo
    import pandas as pd
    from sklearn.cluster import KMeans
    import altair as alt
    return KMeans, alt, mo, pd


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
    return


if __name__ == "__main__":
    app.run()

# marimo_project

## Interactive SMILES Visualization
This project provides an interactive web app to visualize molecular structures derived from a Rhea dataset of SMILES values.

Using dimensionality reduction (UMAP) and clustering (KMeans), the app plots the molecules in 2D space based on structural similarities.
Users can explore the embeddings, view selected 2D structures, and interactively zoom, pan, and select molecules.

## Features
- Interactive scatterplot of UMAP projections of SMILES data.
- Dynamic clustering by molecular features (KMeans).
- On selection, view the 2D molecular structures.
- Smooth zooming, panning, and dynamic table generation.
- Based on SMILES from a Rhea dataset (representing proteins or molecules associated with proteins — biological context under study).

## Built With
Marimo — Python reactive web app framework.

Pandas — data manipulation.

scikit-learn — clustering.

Altair — interactive plotting.

RDKit — molecule generation and visualization.

py3Dmol — 3D molecule rendering (planned extensions).


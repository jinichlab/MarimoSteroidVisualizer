import marimo

__generated_with = "0.10.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np

    #---------------------------------------------
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit import DataStructs
    import rdkit
    return (
        AllChem,
        Chem,
        DataStructs,
        IPythonConsole,
        mo,
        np,
        pd,
        rdFingerprintGenerator,
        rdkit,
    )


@app.cell
def _(mol, pd):
    with open("uniprot/rhea_steroid_smiles.txt", "r") as file:
        rows = []
        for line in file:
            value = line.strip().split(' ')
            rows.append([f"{value[0]} {value[1]}", value[3]])

            #Adding the descriptors (MW, logP) to the descriptors list
            mol(value[3])
    rows

    df = pd.DataFrame(rows, columns=['Compounds', 'SMILES'])
    df.tail()
    return df, file, line, rows, value


@app.cell
def _(rdFingerprintGenerator):
    fmgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2,
                      atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())
    return (fmgen,)


@app.cell
def _(Chem):
    def mol(smile):
        return Chem.MolFromSmiles(smile)
    return (mol,)


@app.cell
def _(df, fmgen, mol):
    cfp = fmgen.GetCountFingerprint(mol(df['SMILES'][0]))
    scfp = fmgen.GetSparseCountFingerprint(mol(df['SMILES'][0]))
    fp = fmgen.GetFingerprint(mol(df['SMILES'][0]))

    #LEARN DIFFERENCE IN THE TYPES OF FINGERPRINTS
    return cfp, fp, scfp


@app.cell
def _(fp):
    fp
    return


@app.cell
def _(cfp, np):
    np.array(cfp)
    return


@app.cell
def _(AllChem, Chem, fmgen, np):
    # Create a MorganGenerator object with specified parameters
    def get_morgan_vector(smiles, radius=2, nBits=1024):
        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to molecule object
        if mol is not None:  # Check if molecule is valid
            # Create the Morgan generator
            morgan_gen = AllChem.rdFingerprintGenerator(radius=radius, nBits=nBits)
            # Generate the Morgan fingerprint
            fingerprint = morgan_gen.GetFingerprint(mol)
            # Convert the fingerprint to numpy array (bit vector)
            return np.array(fmgen(mol))
        else:
            return np.zeros(nBits)
    return (get_morgan_vector,)


@app.cell
def _(AllChem, df, mol, np):
    bit={}
    morganfp2=AllChem.GetMorganFingerprintAsBitVect(mol(df['SMILES'][0]),useChirality=True, radius=2, nBits = 1024, bitInfo=bit)
    mfpvector2 = np.array(morganfp2)
    np.nonzero(mfpvector2)
    return bit, mfpvector2, morganfp2


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

from analysis import MolecularManifold

manifold = MolecularManifold('ras_orbitals.txt')


def plot_orbital_composition(manifold: MolecularManifold):
    '''Plots the orbital composition of each molecular orbital in a manifold'''
    import plotly.graph_objects as go
    import numpy as np
    import plotly.express as px

    df = manifold.to_dataframe(ingore_n=False)
    #plot composition against index ( y-axis = index) as a line plot
    fig = px.line(df, y=df.index, x=df.columns, title='Orbital Composition')
    fig.show()


def plot_orbital_composition_heatmap(manifold: MolecularManifold):
    '''Plots the orbital composition of each molecular orbital in a manifold'''
    import plotly.graph_objects as go
    import plotly.express as px

    df = manifold.to_dataframe(ingore_n=False)

    fig = px.imshow(df, title='Orbital Composition')
    fig.show()

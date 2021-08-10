import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc


adata_selected = False
column_selected = False

#Set options
st.set_option('deprecation.showPyplotGlobalUse', False)

st.title('Cell Visualizer')


uploaded_file = st.file_uploader(label='Upload .h5ad File')

if uploaded_file != None:
    adata = sc.read_h5ad(uploaded_file)
    adata_selected = True
    adata.obs.loc[:,adata.obs.isnull().sum() > 0]

if adata_selected:
    column_selectbox = st.sidebar.selectbox(
        'Column',
        adata.obs.columns
    )

    
    filter_check = st.sidebar.checkbox(label='Filter by category', )

    if filter_check:
        column_values = adata.obs[column_selectbox].tolist()
        column_values = set(column_values)
        filter_selectbox = st.sidebar.selectbox(
            'Filter',
            column_values
    )

        filter_adata = adata[adata.obs[column_selectbox] == filter_selectbox].copy()

        sc.pl.umap(filter_adata, color=[column_selectbox])

    else:
        sc.pl.umap(adata, color=[column_selectbox])

    st.pyplot()



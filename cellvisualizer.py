import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc


#Replace stupid fucking NaN values with 'Unknown'
def replace_nan(column):
    if adata.obs[column].dtype == 'category':
        adata.obs[column] = adata.obs[column].cat.add_categories('Unknown')
        adata.obs[column].fillna('Unknown', inplace =True)
         
def add_filter():
    #Select which column
    column_selectbox = st.sidebar.selectbox(
        'Column ',
        adata.obs.columns
    )
    

    #Replace NaN within column to stop errors
    replace_nan(column_selectbox)

    #Select to filter by category
    filter_check = st.sidebar.checkbox(label='Filter by category ')
    

    #If filtering by category, choose which category
    if filter_check:
        column_values = adata.obs[column_selectbox].tolist()
        column_values = set(column_values)
        filter_selectbox = st.sidebar.selectbox(
            'Filter ',
            column_values
    )
        

        #New adata based on filter
        filter_adata = adata[adata.obs[column_selectbox] == filter_selectbox].copy()
      

        #Filtered umap
        sc.pl.umap(filter_adata, color=[column_selectbox])


    else:
        #Normal umap
        filter_adata = adata
        sc.pl.umap(adata, color=[column_selectbox])

    
adata_selected = False


#Set options
st.set_option('deprecation.showPyplotGlobalUse', False)
st.title('Cell Visualizer')


#Upload .h5ad
uploaded_file = st.file_uploader(label='Upload .h5ad File')


#Check if file was uploaded
if uploaded_file != None:
    original_adata = sc.read_h5ad(uploaded_file)
    adata = original_adata.copy()
    adata_selected = True


if adata_selected:
    add_filter()
    
    #Display umap
    st.pyplot()




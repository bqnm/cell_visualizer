import streamlit as st
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
import scanpy as sc


st.set_page_config(layout="wide")

col1, col2, col3 = st.columns((0.75,2,1.25))


class DataSelection:
    def __init__(self, selection_number, adata):
        self.selection_number = selection_number
        self.adata = adata #.copy()
        self.session_state_index = f'filter_count_{str(selection_number)}'
        if 'filter_count' not in st.session_state:
            st.session_state.filter_count = {}

        if f'filter_count_{str(selection_number)}' not in st.session_state.filter_count:
            st.session_state.filter_count[self.session_state_index] = 0

        self.filter_select_list = adata.obs.columns.tolist() + adata.var_names.tolist()


    def replace_nan(self, column):
        if column.dtype == 'category':
            return column.cat.add_categories('Unknown').fillna('Unknown')
        else:
            return column

    def add_filter(self, filter_count, st_column, key_prefix, adata):
        #Select which column
        column_selectbox = st_column.selectbox(
            f'Column {filter_count}' ,
            self.filter_select_list,
            key=f'{key_prefix} Column {filter_count}'
        )        

        
        #Choose which category / range
        if column_selectbox not in adata.var_names:
            column_series = self.replace_nan(adata.obs[column_selectbox])
            if not is_numeric_dtype(column_series):
                column_values = set(column_series.to_list())
                filter_selectbox = st_column.selectbox(
                    'Filter',
                    column_values,
                    key=f'{key_prefix} Filter {filter_count}'
                )
                    

                #New adata based on filter
                filter_adata = adata[column_series == filter_selectbox] #.copy()
        
                return filter_adata
            
            else:
                min_v, max_v = float(column_series.min()), float(column_series.max())
                filter_range = st_column.slider(
                    label='Filter',
                    min_value=min_v,
                    max_value=max_v,
                    value=(min_v, max_v)
                )
                

                filter_adata = adata[adata.obs[column_selectbox].between(filter_range[0], filter_range[1])]

                return filter_adata

    def create_filters(self):
        col1.title(f'Data Selection {str(self.selection_number)}')
    
        gene_ex_check = col1.checkbox(label='Check for gene expression mode', value=False, key=f'gene_ex_check_{str(self.selection_number)}')

        filter_data = self.adata
        if gene_ex_check == False:
            filter_count = col1.selectbox('Filter count', range(10), key=f'filter_count_box_{str(self.selection_number)}')
            st.session_state.filter_count[self.session_state_index] = filter_count

            for i in range(st.session_state.filter_count[self.session_state_index]):
                filter_data = self.add_filter(i+1, col1, str(self.selection_number), filter_data)
            
            fig = sc.pl.umap(filter_data, show=False, return_fig=True)
            #Display umap
            col2.pyplot(fig)

        else:
            gene_filter_selectbox = col1.selectbox('Gene filter', self.adata.var_names, key=f'gene_filter_selectbox_{self.selection_number}')
            fig = sc.pl.umap(filter_data, show=False, return_fig=True, color=[gene_filter_selectbox])
            col2.pyplot(fig)
            
data_selection_dict = {}

def create_data_selections():
    for i in range(data_selection_count):
        data_selection_dict[f'DataSelection_{i+1}'] = DataSelection(i+1, adata)
        data_selection_dict[f'DataSelection_{i+1}'].create_filters()
    

adata_selected = False


#Set options
col3.title('Cell Visualizer')


#Upload .h5ad
uploaded_file = col3.file_uploader(label='Upload .h5ad File')


#Check if file was uploaded
if uploaded_file != None:
    adata = sc.read_h5ad(uploaded_file)
    adata_selected = True
    
if adata_selected:
    data_selection_count = col1.selectbox('Number of data selections', range(10))
    create_data_selections()
    



    



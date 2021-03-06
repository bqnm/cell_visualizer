import streamlit as st
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
import scanpy as sc
import sys
from pathlib2 import Path

if len(sys.argv) > 1:
    path = Path(sys.argv[1])

else:
    path = Path('./data')

def pick_data_file(container, data_path):
    afiles = list(data_path.glob("**/*.h5ad"))
    num_files = len(afiles)
    options = list(range(len(afiles) + 1))

    def format_path(idx):
        if idx == 0:
            if num_files > 0:
                return "Choose a file..."
            else:
                return f'No files found in "{data_path}"'
        else:
            p = afiles[idx - 1]
            name = p.name
            parent = p.parent.relative_to(data_path).name
            if parent:
                item_str = f"{name} (from {parent})"
            else:
                item_str = name
            return item_str

    item = container.selectbox("File", options, format_func=format_path)

    if item > 0:
        return sc.read_h5ad(afiles[item - 1])
    else:
        return None

st.set_page_config(layout="wide")               
col1, col2 = st.columns((0.75,2))


class DataSelection:
    def __init__(self, selection_number, adata):
        self.selection_number = selection_number
        self.adata = adata
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

    def add_filter(self, filter_count, st_column, key_prefix, adata, filter_type):
        column_selectbox = st_column.selectbox(
            f'Column {filter_count}' ,
            self.filter_select_list,
            key=f'{key_prefix} Column {filter_count}'
            )
        if column_selectbox not in adata.var_names:
            column_series = adata.obs[column_selectbox]
            if not is_numeric_dtype(column_series):
                column_series = self.replace_nan(adata.obs[column_selectbox])
                column_values = set(column_series.to_list())
                filter_selectbox = st_column.selectbox(
                    'Filter',
                    column_values,
                    key=f'{key_prefix} Filter {filter_count}'
                )
                
                filter_series = column_series == filter_selectbox
                filter_label = column_selectbox
                return filter_series
            else:
                min_v, max_v = float(column_series.min()), float(column_series.max())
                filter_range = st_column.slider(
                    label='Filter',
                    min_value=min_v,
                    max_value=max_v,
                    value=(min_v, max_v),
                    key=f'{key_prefix} Slider {filter_count}'
                )
                
                filter_series = adata.obs[column_selectbox].between(filter_range[0], filter_range[1])
                return filter_series

        else:
            gene_vals = adata[:,column_selectbox].X
            min_v, max_v = float(gene_vals.min()), float(gene_vals.max())
            filter_range = st_column.slider(
                    label='Filter',
                    min_value=min_v,
                    max_value=max_v,
                    value=(min_v, max_v),
                    key=f'{key_prefix} Slider {filter_count}'
                )

            filter_series = (gene_vals > filter_range[0]) & (gene_vals < filter_range[1])
            return filter_series

    def add_color_filter(self, st_column, key_prefix):
        column_selectbox = st_column.selectbox(
            'Color column' ,
            self.filter_select_list,
            key=f'{key_prefix} Color Column'
            )

        return column_selectbox

    def create_filters(self):
        col1.title(f'Data Selection {str(self.selection_number)}')

        filter_adata = self.adata
        
        filter_count = col1.number_input('Filter count', min_value=0, key=f'filter_count_box_{str(self.selection_number)}')
        st.session_state.filter_count[self.session_state_index] = filter_count

        for i in range(int(st.session_state.filter_count[self.session_state_index])):
            filter_series = self.add_filter(i+1, col1, str(self.selection_number), filter_adata, 'selection')
            filter_adata = filter_adata[filter_series]
        
        color_checkbox = col1.checkbox('Color map', key=f'color_checkbox_{str(self.selection_number)}')

        if color_checkbox:
            color_filter = self.add_color_filter(col1, str(self.selection_number))
            fig = sc.pl.umap(filter_adata, show=False, return_fig=True, color=[color_filter])
        
        else:
            fig = sc.pl.umap(filter_adata, show=False, return_fig=True)
        
        col2.pyplot(fig)
            
data_selection_dict = {}

def create_data_selections():
    for i in range(data_selection_count):
        data_selection_dict[f'DataSelection_{i+1}'] = DataSelection(i+1, adata)
        data_selection_dict[f'DataSelection_{i+1}'].create_filters()



col1.title('Cell Visualizer')

default_dataset = col1.checkbox('Use pbmc3k_processed')

if default_dataset:
    adata = sc.datasets.pbmc3k_processed()
else:  
    adata = pick_data_file(col1, path)
    
if adata:
    data_selection_count = col1.number_input('Number of data selections', min_value=0)
    create_data_selections()
    



    



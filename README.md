# Cell Visualizer

For looking at and plotting .h5ad files created in [Scanpy](https://scanpy.readthedocs.io/en/stable/)

You can access the visualizer online [here](https://share.streamlit.io/bqnm/cell_visualizer/cellvisualizer.py), however if you want to use your own .h5ad file other than the default data provided (pbmc3k_processed) you need to run it locally.

To use the visualizer locally, while in the cell_visualizer directory run: 

'streamlit run cellvisualizer.py --server.maxUploadSize=(MAX UPLOAD SIZE REQUIRED HERE)' 


Unless you provide a different folder path as a command line argument this will automatically search for a folder called 'data'. This is where you should put your .h5ead files.

https://user-images.githubusercontent.com/365396/130834912-077ab0c1-7cca-460a-9476-cc9a85436da2.mp4



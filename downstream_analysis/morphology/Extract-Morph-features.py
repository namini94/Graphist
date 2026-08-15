# import module
import stlearn as st
from pathlib import Path
import pandas as pd
#st.settings.set_figure_params(dpi=180)


# specify PATH to data
BASE_PATH = Path("/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/outs")

# spot tile is the intermediate result of image pre-processing
TILE_PATH = Path("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/tmp/tiles")
TILE_PATH.mkdir(parents=True, exist_ok=True)



# load data
data = st.Read10X(BASE_PATH)


# pre-processing for gene count table
st.pp.filter_genes(data,min_cells=1)
st.pp.normalize_total(data)
st.pp.log1p(data)


# pre-processing for spot image
st.pp.tiling(data, TILE_PATH)

# this step uses deep learning model to extract high-level features from tile images
# may need few minutes to be completed

import ssl
ssl._create_default_https_context = ssl._create_unverified_context
st.pp.extract_feature(data, cnn_base = "xception")
#st.pp.extract_feature(data, cnn_base = "vgg16")
#st.pp.extract_feature(data, cnn_base = "resnet50")

morph = data.obsm['X_morphology']
morph = pd.DataFrame(morph)
morph.index = data.obs_names

pd.DataFrame(morph).to_csv("/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/Morph-XCEPTION.csv",index=True)

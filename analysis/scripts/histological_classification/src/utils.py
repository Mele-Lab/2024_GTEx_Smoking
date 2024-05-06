import os
from histolab.slide import Slide
from pathlib import Path
import pandas as pd
import numpy as np

import yaml
def get_parameters():
    with open("./conf/parameters.yml") as params:
        params_dict = yaml.safe_load(params)
    return params_dict


params = get_parameters()

path__gtex_portal_data = params["PATH__GTExPortalDB"]
path__smoker_annotation = params["PATH__SmokerStatus"]

dfpd_GtexPortal = pd.read_csv(path__gtex_portal_data).rename(
        {
            "Tissue Sample ID":"sample_id", 
            "Subject ID":"subject_id"
        }, axis=1
    )

dfpd_SmokerStatus = pd.read_csv(path__smoker_annotation, sep=";").rename(
        {
            "Donor":"subject_id", 
            "SmokerStatus":"smoker_status"
        }, axis=1
    )


def get_pathology_array(sid, tissue_type, dfpd_GtexPortal=dfpd_GtexPortal):
    dfpd_pathologies = dfpd_GtexPortal.query(f"Tissue == '{tissue_type}'")["Pathology Categories"].str.split(", ").explode()
    dfpd_pathologies = pd.get_dummies(dfpd_pathologies).groupby(level=0).sum()
    dfpd_pathologies = dfpd_pathologies.set_index(dfpd_GtexPortal.loc[dfpd_pathologies.index]["sample_id"])
    return dfpd_pathologies.loc[sid].to_dict()

def get_smoker_status(subject_id):
    dfpd_temp = dfpd_SmokerStatus[dfpd_SmokerStatus["subject_id"] == subject_id]
    if dfpd_temp.empty:
        return "unknown"
    elif dfpd_temp.MHSMKPRD.to_list()[0] in 'Non Smoker':
        return "Non Smoker"
    elif dfpd_temp.MHSMKCMT.to_list()[0] in 'Non Smoker':
        return "Smoker"
    else:
        smoker_status = dfpd_temp.smoker_status.fillna("unknown").to_list()[0]
        return smoker_status

def get_age(subject_id):
    dfpd_temp = dfpd_SmokerStatus[dfpd_SmokerStatus["subject_id"] == subject_id]
    age = dfpd_temp["age"].to_list()[0]
    return age

class GtexSample():
    def __init__(self, path__gtex_portal_data=path__gtex_portal_data):
        self.database = pd.read_csv(path__gtex_portal_data).rename(
                {
                    "Tissue Sample ID":"sample_id", 
                    "Subject ID":"subject_id"
                }, axis=1
            )
        
        self.attributes = {
            "tissues": list(self.database["Tissue"].unique())
        }

    def get_subsample(self, key, value):
        df_temp = self.database[self.database[key] == value]
        lt_sid = list(df_temp.sample_id.unique())

        print(f"[+] Found {len(lt_sid)} subjects for {value}")
        return list(df_temp.sample_id.unique())

    def compare_subsamples(self, key, values, return_subsamples=False):
        dct_sets = {t:set(["-".join(i.split("-")[:-1]) for i in self.get_subsample(key, t)]) for t in values}

        import matplotlib.pyplot as plt
        import venn

        venn_diagram = venn.venn(dct_sets, cmap="plasma", fontsize=8, figsize=(5,5))
        plt.title(f'Number of subjects per {key}')
        plt.show()

        if return_subsamples:
            return dct_sets


class GtexSlide():
    def __init__(self, sid, output):
        self.sid = sid
        self.sample_id = "-".join(sid.split('-')[:-1])
        self.output = output
        self.tile_path = f'{self.output}/{self.sid}'
        self.slide_path = f'{self.output}/{self.sid}.svs'
        self.attributes = {
            "GTExPortal": dfpd_GtexPortal[dfpd_GtexPortal["sample_id"] == self.sid].to_dict("records")[0]
        }

        self.tissue = self.attributes["GTExPortal"]["Tissue"]

        self.attributes.update(
            {"pathologies" : get_pathology_array(self.sid, self.tissue)}
            )

        self.attributes.update(
            {"smoker_status": get_smoker_status(self.sample_id)}
        )

        if not os.path.exists(self.tile_path):
            Path(self.tile_path).mkdir(parents=True, exist_ok=True)
            print("[+] Tile path created.")
        if not os.path.exists(self.output):
            Path(self.output).mkdir(parents=True, exist_ok=True)
            print("[+] Output path created.")

        self.slide = Slide(self.slide_path, processed_path=self.tile_path, use_largeimage=False)

    def delete_wsi(self):
        rm_svs = f"rm {self.slide_path}"
        os.system(rm_svs)
        print(f"{self.sid} removed")
        print()
    
    def get_wsi(self):
        # Check if the file already exists
        if os.path.exists(self.slide_path):
            print(f" [-] File {self.sid} already exists.")
            return False
        else:
            try:
                print(f"[+] Downloading wsi {self.sid}...")
                curl_cmd = f"curl 'https://brd.nci.nih.gov/brd/imagedownload/{self.sid}' --output '{self.slide_path}'"
                os.system(curl_cmd)
                return True
            except Exception as e:
                print(f'Problems downloading wsi {self.sid}: {e}')
    
    def info(self):
        print(f"Slide name: {self.slide.name}")
        print(f"Levels: {self.slide.levels}")
        print("Dimensions and magnification at:")
        print(f"  - level 0: {self.slide.dimensions}, {self.slide.level_magnification_factor(level=0)}")
        print(f"  - level 1: {self.slide.level_dimensions(level=1)}, {self.slide.level_magnification_factor(level=1)}")
        print(f"  - level 2: {self.slide.level_dimensions(level=2)}, {self.slide.level_magnification_factor(level=2)}")
        print(
        "Native magnification factor:",
        self.slide.level_magnification_factor()
        )

    def break_into_tiles(self, level=1, size=(512,512), tissue_percent=85, thumb=False):
        from histolab.masks import BiggestTissueBoxMask, TissueMask
        from histolab.filters.image_filters import OtsuThreshold
        #bigbox_mask = BiggestTissueBoxMask(OtsuThreshold())
        tissue_mask = TissueMask(OtsuThreshold())

        mask = tissue_mask


        # Get tiles: With NucleiScorer
        from histolab.tiler import ScoreTiler, GridTiler
        from histolab.scorer import NucleiScorer, CellularityScorer

        self.prefix_path = f"sz{size[0]}_lv{level}_tp{tissue_percent}/tiles__raw/"

        tiles_extractor = GridTiler(
                #scorer = NucleiScorer(),
                tile_size=size,
                #n_tiles=20,
                level=level,
                check_tissue=True,
                tissue_percent=tissue_percent,
                pixel_overlap=0,
                #mpp=494.2,
                prefix=self.prefix_path,
                suffix=".png"
            )
        
        tiles_extractor.extract(self.slide, extraction_mask=mask)#, report_path=SUMMARY_PATH)

        if thumb:
            img = tiles_extractor.locate_tiles(slide=self.slide, extraction_mask=mask, outline="black", linewidth=1)
            img.save(f"{self.tile_path}/sz{size[0]}_lv{level}_tp{tissue_percent}/thumb_extraction.png")
        
    def get_tiles_path(self, prefix_path=None):
        if prefix_path:
            self.prefix_path = f"{self.tile_path}/{prefix_path}"
        
        path = f"{self.prefix_path}"
        return os.listdir(path)
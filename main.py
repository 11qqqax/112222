import argparse
import yaml

from lib.runner import Runner

# cd D:\mesh_data\qh\IntrinsicGarmAlign-main
# python main.py configs\g01.yaml

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("cfgfile", type=str)
    # args = parser.parse_args()

    # 直接指定配置文件路径
    cfgfile_path = "configs\g01.yaml"

    with open(cfgfile_path, "r") as cfgfile:
        cfg = yaml.safe_load(cfgfile)

    runner = Runner(cfg)
    runner.compute_eigenfuncs()
    runner.set_template()
    runner.pose_by_lbs()
    runner.deform_coarse_stage()
    runner.rectify_embeddings_by_FM()
    runner.deform_refining_stage()
    runner.shape_transfer()

Rscript src/check_dependencies.R

# load conda env and set path to it
conda env create -f env.yml
conda info --envs | grep -Po "dipdsa\K.*" | sed 's: ::g' > conda_path.txt

